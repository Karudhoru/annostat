"""Deterministic structural validation for paired prokaryotic FASTA/GFF3 files.

Validation is deliberately limited to properties that can be established from
the two input files.  Biological expectations that can have legitimate
exceptions belong in :mod:`annostat.qc` and are reported as inspection findings.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from pathlib import Path

from annostat import __version__
from annostat.models import Feature
from annostat.parsers import parse_fasta, parse_gff


VALIDATION_SCHEMA_VERSION = "1.0"
RULESET_VERSION = "2026.08.2"
_IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")


@dataclass(frozen=True, slots=True)
class ValidationRule:
    """Describe one stable validation rule and its scientific provenance."""

    rule_id: str
    title: str
    category: str
    default_severity: str
    rationale: str
    source_url: str


@dataclass(frozen=True, slots=True)
class ValidationFinding:
    """Represent one deterministic, machine-readable validation finding."""

    rule_id: str
    severity: str
    message: str
    seqid: str = ""
    feature_id: str = ""
    start: int | None = None
    end: int | None = None
    observed: str = ""
    expected: str = ""

    def as_dict(self) -> dict[str, str | int | None]:
        """Return a stable JSON/TSV representation."""

        return asdict(self)


@dataclass(frozen=True, slots=True)
class ValidatedAnnotation:
    """Hold parsed inputs together with their deterministic validation result.

    Analysis workflows use this internal bundle to avoid parsing the same FASTA
    and GFF3 files again after validation. ``validate_annotation`` remains the
    stable JSON-safe public interface for callers that only need validation.
    """

    genome: dict[str, str]
    features: tuple[Feature, ...]
    validation: dict[str, object]


RULES = {
    rule.rule_id: rule
    for rule in (
        ValidationRule(
            "GFF_VERSION",
            "GFF3 version directive",
            "format",
            "warning",
            "A GFF3 document should declare version 3 near the beginning of the file.",
            "https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md",
        ),
        ValidationRule(
            "GFF_PARSE",
            "Parseable GFF3 records",
            "format",
            "error",
            "Every feature record must contain valid values in the nine GFF3 columns.",
            "https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md",
        ),
        ValidationRule(
            "FASTA_PARSE",
            "Parseable FASTA records",
            "format",
            "error",
            "Each FASTA sequence must have a non-empty, unique identifier and sequence.",
            "https://www.ncbi.nlm.nih.gov/genbank/fastaformat/",
        ),
        ValidationRule(
            "FASTA_ALPHABET",
            "IUPAC DNA alphabet",
            "sequence",
            "error",
            "Genome sequence characters must be valid IUPAC DNA symbols.",
            "https://www.ncbi.nlm.nih.gov/IEB/ToolBox/SDKDOCS/BIOSEQ.HTML",
        ),
        ValidationRule(
            "SEQID_MATCH",
            "GFF3 and FASTA sequence identifiers agree",
            "cross_file",
            "error",
            "A feature can only be interpreted against a FASTA record with the same identifier.",
            "https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md",
        ),
        ValidationRule(
            "COORDINATE_BOUNDS",
            "Feature coordinates fit the reference",
            "cross_file",
            "error",
            "Linear feature coordinates must lie within the referenced sequence.",
            "https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md",
        ),
        ValidationRule(
            "PARENT_EXISTS",
            "Parent identifiers resolve",
            "relationships",
            "error",
            "Every Parent value must identify a feature in the same GFF3 document.",
            "https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md",
        ),
        ValidationRule(
            "ID_CONSISTENCY",
            "Repeated IDs describe one multipart feature",
            "relationships",
            "error",
            "Repeated IDs are valid only for rows representing parts of the same feature.",
            "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/",
        ),
        ValidationRule(
            "CDS_PHASE",
            "CDS phases preserve the reading frame",
            "translation",
            "error",
            "CDS phase specifies how many bases precede the next complete codon in biological order.",
            "https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md",
        ),
        ValidationRule(
            "TRANSLATION_TABLE",
            "Annotation translation table is supported and consistent",
            "translation",
            "error",
            "NCBI GFF3 CDS or region rows may declare transl_table; one consistent supported table is required for deterministic translation.",
            "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/",
        ),
    )
}


def _sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest for provenance."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _finding(
    rule_id: str,
    message: str,
    *,
    feature: Feature | None = None,
    severity: str | None = None,
    observed: object = "",
    expected: object = "",
) -> ValidationFinding:
    """Construct a finding using the rule's default severity."""

    return ValidationFinding(
        rule_id=rule_id,
        severity=severity or RULES[rule_id].default_severity,
        message=message,
        seqid=feature.seqid if feature else "",
        feature_id=feature.id if feature else "",
        start=feature.start if feature else None,
        end=feature.end if feature else None,
        observed=str(observed),
        expected=str(expected),
    )


def _version_findings(path: Path) -> list[ValidationFinding]:
    """Check the GFF version directive without requiring successful parsing."""

    declared = ""
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("##gff-version"):
                declared = line.removeprefix("##gff-version").strip()
            break
    if declared == "3":
        return []
    if declared:
        return [_finding("GFF_VERSION", "unsupported GFF version", observed=declared, expected="3")]
    return [_finding("GFF_VERSION", "missing ##gff-version 3 directive", expected="3")]


def _circular_sequence_ids(features: Iterable[Feature]) -> frozenset[str]:
    """Return sequence IDs explicitly marked circular by a region feature."""

    return frozenset(
        feature.seqid
        for feature in features
        if feature.type == "region"
        and feature.attributes.get("Is_circular", "").strip().lower() == "true"
    )


def _alphabet_findings(genome: Mapping[str, str]) -> list[ValidationFinding]:
    """Report non-IUPAC characters without treating ambiguity codes as errors."""

    findings = []
    for seqid, sequence in sorted(genome.items()):
        invalid = sorted(set(sequence.upper()) - _IUPAC_DNA)
        if invalid:
            findings.append(
                ValidationFinding(
                    rule_id="FASTA_ALPHABET",
                    severity="error",
                    message="sequence contains non-IUPAC DNA characters",
                    seqid=seqid,
                    observed="".join(invalid),
                    expected="ACGTRYSWKMBDHVN",
                )
            )
    return findings


def _identifier_findings(features: list[Feature]) -> list[ValidationFinding]:
    """Validate Parent references and the consistency of repeated feature IDs."""

    by_id: dict[str, list[Feature]] = defaultdict(list)
    for feature in features:
        explicit_id = feature.attributes.get("ID", "").strip()
        if explicit_id:
            by_id[explicit_id].append(feature)

    findings: list[ValidationFinding] = []
    identifiers = set(by_id)
    for feature in features:
        for parent in filter(None, feature.attributes.get("Parent", "").split(",")):
            if parent not in identifiers:
                findings.append(
                    _finding(
                        "PARENT_EXISTS",
                        f"Parent {parent!r} does not resolve to a feature ID",
                        feature=feature,
                        observed=parent,
                    )
                )

    for identifier, parts in sorted(by_id.items()):
        if len(parts) < 2:
            continue
        signatures = {(part.seqid, part.type, part.strand) for part in parts}
        parents = {part.attributes.get("Parent", "") for part in parts}
        if len(signatures) != 1 or len(parents) != 1:
            findings.append(
                _finding(
                    "ID_CONSISTENCY",
                    "repeated ID has inconsistent seqid, type, strand, or Parent",
                    feature=parts[0],
                    observed=identifier,
                    expected="identical seqid, type, strand, and Parent",
                )
            )
    return findings


def _coordinate_findings(
    features: list[Feature],
    genome: Mapping[str, str],
    circular_seqids: frozenset[str],
) -> list[ValidationFinding]:
    """Validate cross-file references and coordinate bounds."""

    findings = []
    for feature in features:
        if feature.seqid not in genome:
            findings.append(
                _finding(
                    "SEQID_MATCH",
                    "feature sequence ID is absent from FASTA",
                    feature=feature,
                    observed=feature.seqid,
                )
            )
            continue
        sequence_length = len(genome[feature.seqid])
        if feature.seqid in circular_seqids and feature.start > sequence_length:
            findings.append(
                _finding(
                    "COORDINATE_BOUNDS",
                    "circular feature starts beyond the real sequence coordinate range",
                    feature=feature,
                    observed=feature.start,
                    expected=f"<= {sequence_length}",
                )
            )
        if feature.end > sequence_length and feature.seqid not in circular_seqids:
            findings.append(
                _finding(
                    "COORDINATE_BOUNDS",
                    "feature extends beyond a linear FASTA sequence",
                    feature=feature,
                    observed=feature.end,
                    expected=f"<= {sequence_length}",
                )
            )
        if feature.seqid in circular_seqids and feature.length > sequence_length:
            findings.append(
                _finding(
                    "COORDINATE_BOUNDS",
                    "circular feature spans more than one complete sequence length",
                    feature=feature,
                    observed=feature.length,
                    expected=f"<= {sequence_length}",
                )
            )
    return findings


def _phase_findings(features: list[Feature]) -> list[ValidationFinding]:
    """Validate phase continuity for single- and multipart CDS features."""

    groups: dict[str, list[Feature]] = defaultdict(list)
    for index, feature in enumerate(features):
        if feature.type != "CDS":
            continue
        identifier = feature.attributes.get("ID") or f"__row_{index}"
        groups[identifier].append(feature)

    findings = []
    for identifier, parts in sorted(groups.items()):
        ordered = sorted(parts, key=lambda item: item.start, reverse=parts[0].strand == "-")
        first_phase = ordered[0].phase
        if first_phase is None:
            findings.append(
                _finding(
                    "CDS_PHASE",
                    "CDS phase is missing",
                    feature=ordered[0],
                    observed=".",
                    expected="0, 1, or 2",
                )
            )
            continue
        coding_bases = ordered[0].length - first_phase
        for part in ordered[1:]:
            expected_phase = (3 - coding_bases % 3) % 3
            if part.phase != expected_phase:
                findings.append(
                    _finding(
                        "CDS_PHASE",
                        f"multipart CDS {identifier!r} does not preserve its reading frame",
                        feature=part,
                        observed="." if part.phase is None else part.phase,
                        expected=expected_phase,
                    )
                )
            coding_bases += part.length
    return findings


def _translation_table_findings(features: list[Feature]) -> list[ValidationFinding]:
    """Report malformed, unsupported, or conflicting CDS translation tables."""

    findings: list[ValidationFinding] = []
    declared: dict[int, Feature] = {}
    expected = "4, 11, or 25"
    for feature in features:
        if feature.type not in {"CDS", "region"} or "transl_table" not in feature.attributes:
            continue
        raw_code = feature.attributes["transl_table"].strip()
        try:
            code = int(raw_code)
        except ValueError:
            findings.append(
                _finding(
                    "TRANSLATION_TABLE",
                    "feature transl_table is not an integer",
                    feature=feature,
                    observed=raw_code,
                    expected=expected,
                )
            )
            continue
        if code not in {4, 11, 25}:
            findings.append(
                _finding(
                    "TRANSLATION_TABLE",
                    "feature declares an unsupported translation table",
                    feature=feature,
                    observed=code,
                    expected=expected,
                )
            )
            continue
        declared.setdefault(code, feature)
    if len(declared) > 1:
        choices = ", ".join(map(str, sorted(declared)))
        findings.append(
            _finding(
                "TRANSLATION_TABLE",
                "annotation features declare conflicting translation tables",
                feature=declared[min(declared)],
                observed=choices,
                expected="one consistent table",
            )
        )
    return findings


def _sort_key(finding: ValidationFinding) -> tuple[object, ...]:
    """Return the canonical ordering used by every output format."""

    severity_order = {"error": 0, "warning": 1, "info": 2}
    return (
        severity_order.get(finding.severity, 9),
        finding.rule_id,
        finding.seqid,
        finding.start if finding.start is not None else -1,
        finding.end if finding.end is not None else -1,
        finding.feature_id,
        finding.message,
    )


def load_and_validate_annotation(
    fasta_path: Path,
    gff_path: Path,
) -> ValidatedAnnotation:
    """Parse paired inputs once and return their validation and loaded data.

    ``valid`` reflects structural errors only. Warnings identify portability or
    metadata concerns and do not assert that the underlying biology is wrong.
    """

    fasta_path = Path(fasta_path)
    gff_path = Path(gff_path)
    findings: list[ValidationFinding] = []
    genome: dict[str, str] = {}
    features: list[Feature] = []

    try:
        genome = parse_fasta(fasta_path)
    except (OSError, ValueError) as error:
        findings.append(_finding("FASTA_PARSE", str(error)))
    try:
        findings.extend(_version_findings(gff_path))
        features = list(parse_gff(gff_path))
    except (OSError, ValueError) as error:
        findings.append(_finding("GFF_PARSE", str(error)))

    if genome:
        findings.extend(_alphabet_findings(genome))
    if features:
        circular_seqids = _circular_sequence_ids(features)
        findings.extend(_identifier_findings(features))
        if genome:
            findings.extend(_coordinate_findings(features, genome, circular_seqids))
        findings.extend(_phase_findings(features))
        findings.extend(_translation_table_findings(features))

    findings.sort(key=_sort_key)
    severity_counts = Counter(finding.severity for finding in findings)
    rule_counts = Counter(finding.rule_id for finding in findings)
    input_files = {"fasta": str(fasta_path), "gff3": str(gff_path)}
    input_sha256 = {}
    for name, path in (("fasta", fasta_path), ("gff3", gff_path)):
        if path.is_file():
            input_sha256[name] = _sha256(path)

    result = {
        "schema_version": VALIDATION_SCHEMA_VERSION,
        "ruleset_version": RULESET_VERSION,
        "annostat_version": __version__,
        "valid": severity_counts["error"] == 0,
        "input_files": input_files,
        "input_sha256": input_sha256,
        "record_counts": {
            "fasta_sequences": len(genome),
            "gff3_rows": len(features),
            "cds_rows": sum(feature.type == "CDS" for feature in features),
        },
        "severity_counts": {
            severity: severity_counts.get(severity, 0)
            for severity in ("error", "warning", "info")
        },
        "rule_counts": dict(sorted(rule_counts.items())),
        "findings": [finding.as_dict() for finding in findings],
        "rules": [asdict(RULES[rule_id]) for rule_id in sorted(RULES)],
    }
    fingerprint_payload = {
        key: result[key]
        for key in (
            "schema_version", "ruleset_version", "input_sha256", "record_counts",
            "severity_counts", "rule_counts", "findings",
        )
    }
    result["scientific_fingerprint"] = hashlib.sha256(
        json.dumps(
            fingerprint_payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode("utf-8")
    ).hexdigest()
    return ValidatedAnnotation(
        genome=genome,
        features=tuple(features),
        validation=result,
    )


def validate_annotation(fasta_path: Path, gff_path: Path) -> dict[str, object]:
    """Validate paired inputs and return a deterministic, JSON-safe result."""

    return load_and_validate_annotation(fasta_path, gff_path).validation


def canonical_json(result: Mapping[str, object]) -> str:
    """Serialize a validation result reproducibly across repeated runs."""

    return json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def write_validation(output_dir: Path, result: Mapping[str, object]) -> list[Path]:
    """Write canonical JSON plus a flat TSV finding table."""

    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "validation.json"
    table_path = output_dir / "validation.tsv"
    json_path.write_text(canonical_json(result), encoding="utf-8")
    fieldnames = [
        "rule_id", "severity", "message", "seqid", "feature_id", "start",
        "end", "observed", "expected",
    ]
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(result["findings"])
    return [json_path, table_path]
