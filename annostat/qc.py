"""Internal consistency and annotation-quality checks for bacterial genomes."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass

from annostat.models import CdsSequence, Feature


_STANDARD_RRNA_TYPES = ("5S", "16S", "23S")
_STANDARD_TRNA_AMINO_ACIDS = frozenset(
    {
        "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
        "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
    }
)
_TRNA_PATTERN = re.compile(r"tRNA-([A-Za-z]{3})", re.IGNORECASE)
_START_CODONS = {
    4: frozenset({"TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}),
    11: frozenset({"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}),
    25: frozenset({"TTG", "ATG", "GTG"}),
}
_STOP_CODONS = {
    4: frozenset({"TAA", "TAG"}),
    11: frozenset({"TAA", "TAG", "TGA"}),
    25: frozenset({"TAA", "TAG"}),
}


@dataclass(frozen=True, slots=True)
class AnnotationIssue:
    """Represent one reproducible annotation-quality finding."""

    issue_type: str
    severity: str
    seqid: str = ""
    feature_id: str = ""
    related_feature_id: str = ""
    start: int | None = None
    end: int | None = None
    strand: str = ""
    details: str = ""

    def as_dict(self) -> dict[str, str | int | None]:
        """Return the finding as a dictionary suitable for CSV output."""

        return asdict(self)


def _attribute_is_true(feature: Feature, name: str) -> bool:
    """Interpret common truthy spellings in a GFF3 attribute."""

    return feature.attributes.get(name, "").strip().lower() in {"1", "true", "yes"}


def _is_pseudogene(feature: Feature) -> bool:
    """Return whether GFF3 type or attributes explicitly mark a pseudogene."""

    return (
        "pseudogene" in feature.type.lower()
        or "pseudogene" in feature.attributes
        or _attribute_is_true(feature, "pseudo")
    )


def _is_partial(feature: Feature) -> bool:
    """Return whether a feature carries an explicit partial-boundary marker."""

    return (
        _attribute_is_true(feature, "partial")
        or "start_range" in feature.attributes
        or "end_range" in feature.attributes
    )


def _is_five_prime_partial(feature: Feature) -> bool:
    """Return whether the biological CDS start is explicitly incomplete."""

    boundary = "end_range" if feature.strand == "-" else "start_range"
    if boundary in feature.attributes:
        return True
    return _attribute_is_true(feature, "partial") and not (
        {"start_range", "end_range"} & feature.attributes.keys()
    )


def _is_three_prime_partial(feature: Feature) -> bool:
    """Return whether the biological CDS end is explicitly incomplete."""

    boundary = "start_range" if feature.strand == "-" else "end_range"
    if boundary in feature.attributes:
        return True
    return _attribute_is_true(feature, "partial") and not (
        {"start_range", "end_range"} & feature.attributes.keys()
    )


def _finding_for_feature(
    issue_type: str,
    severity: str,
    feature: Feature,
    details: str,
    related_feature_id: str = "",
) -> AnnotationIssue:
    """Build a quality finding populated with one feature's coordinates."""

    return AnnotationIssue(
        issue_type=issue_type,
        severity=severity,
        seqid=feature.seqid,
        feature_id=feature.id,
        related_feature_id=related_feature_id,
        start=feature.start,
        end=feature.end,
        strand=feature.strand,
        details=details,
    )


def _circular_intervals(feature: Feature, sequence_length: int) -> list[tuple[int, int]]:
    """Represent one circular feature as one or two normalized intervals."""

    if feature.length >= sequence_length:
        return [(1, sequence_length)]
    start = (feature.start - 1) % sequence_length + 1
    end = start + feature.length - 1
    if end <= sequence_length:
        return [(start, end)]
    return [(start, sequence_length), (1, end - sequence_length)]


def _overlap_findings(
    cds_features: Iterable[Feature],
    genome: Mapping[str, str] | None = None,
    circular_seqids: frozenset[str] = frozenset(),
) -> list[AnnotationIssue]:
    """Find overlapping and fully contained CDS pairs with a sweep-line pass."""

    by_sequence: dict[str, list[Feature]] = defaultdict(list)
    for feature in cds_features:
        by_sequence[feature.seqid].append(feature)

    findings: list[AnnotationIssue] = []
    for seqid, sequence_features in by_sequence.items():
        if genome is not None and seqid in circular_seqids and seqid in genome:
            sequence_length = len(genome[seqid])
            ordered = sorted(sequence_features, key=lambda item: (item.start, item.end, item.id))
            segments = sorted(
                (start, end, index)
                for index, feature in enumerate(ordered)
                for start, end in _circular_intervals(feature, sequence_length)
            )
            active_segments: list[tuple[int, int, int]] = []
            overlap_by_pair: Counter[tuple[int, int]] = Counter()
            for start, end, index in segments:
                active_segments = [segment for segment in active_segments if segment[1] >= start]
                for _, other_end, other_index in active_segments:
                    if index == other_index or ordered[index].id == ordered[other_index].id:
                        continue
                    pair = tuple(sorted((index, other_index)))
                    overlap_by_pair[pair] += min(end, other_end) - start + 1
                active_segments.append((start, end, index))
            for (left_index, right_index), overlap_bases in sorted(overlap_by_pair.items()):
                other, current = ordered[left_index], ordered[right_index]
                contained = overlap_bases >= min(
                    current.length, other.length, sequence_length
                )
                relation = "same strand" if other.strand == current.strand else "opposite strands"
                findings.append(
                    _finding_for_feature(
                        "contained_cds" if contained else "cds_overlap",
                        "warning" if contained else "info",
                        current,
                        f"{overlap_bases} bp overlap on {relation}",
                        related_feature_id=other.id,
                    )
                )
            continue
        active: list[Feature] = []
        for current in sorted(sequence_features, key=lambda item: (item.start, item.end, item.id)):
            active = [other for other in active if other.end >= current.start]
            for other in active:
                # Multi-part CDS rows may share an ID and are one biological feature.
                if other.id == current.id:
                    continue
                overlap_bases = min(other.end, current.end) - current.start + 1
                current_contained = current.end <= other.end
                other_contained = (
                    other.start >= current.start and other.end <= current.end
                )
                contained = current_contained or other_contained
                issue_type = "contained_cds" if contained else "cds_overlap"
                severity = "warning" if contained else "info"
                relation = "same strand" if other.strand == current.strand else "opposite strands"
                findings.append(
                    _finding_for_feature(
                        issue_type,
                        severity,
                        current,
                        f"{overlap_bases} bp overlap on {relation}",
                        related_feature_id=other.id,
                    )
                )
            active.append(current)
    return findings


def _annotation_label(feature: Feature) -> tuple[str, str] | None:
    """Return a conservative gene or product label for adjacency comparison."""

    gene = feature.attributes.get("gene", "").strip()
    if gene:
        return "gene", gene.casefold()
    product = feature.attributes.get("product", "").strip()
    if product and product.casefold() != "hypothetical protein":
        return "product", product.casefold()
    return None


def _adjacent_duplicate_findings(cds_features: Iterable[Feature]) -> list[AnnotationIssue]:
    """Report neighboring CDS entries with the same specific annotation label."""

    by_sequence: dict[str, list[Feature]] = defaultdict(list)
    for feature in cds_features:
        by_sequence[feature.seqid].append(feature)

    findings: list[AnnotationIssue] = []
    for sequence_features in by_sequence.values():
        ordered = sorted(sequence_features, key=lambda item: (item.start, item.end, item.id))
        for left, right in zip(ordered, ordered[1:]):
            left_label = _annotation_label(left)
            right_label = _annotation_label(right)
            if left_label is None or left_label != right_label:
                continue
            gap = right.start - left.end - 1
            findings.append(
                _finding_for_feature(
                    "adjacent_duplicate_annotation",
                    "info",
                    right,
                    f"same {right_label[0]} annotation; inter-feature gap {gap} bp",
                    related_feature_id=left.id,
                )
            )
    return findings


def _structural_rna_findings(features: Iterable[Feature]) -> list[AnnotationIssue]:
    """Check for the standard bacterial rRNAs and tRNA amino-acid coverage."""

    observed_rrna: set[str] = set()
    observed_trna: set[str] = set()
    for feature in features:
        annotation = " ".join(
            (
                feature.attributes.get("Name", ""),
                feature.attributes.get("product", ""),
            )
        )
        if feature.type == "rRNA":
            observed_rrna.update(kind for kind in _STANDARD_RRNA_TYPES if kind in annotation)
        if feature.type == "tRNA":
            match = _TRNA_PATTERN.search(annotation)
            if match:
                amino_acid = match.group(1).title()
                if amino_acid in _STANDARD_TRNA_AMINO_ACIDS:
                    observed_trna.add(amino_acid)

    findings = [
        AnnotationIssue(
            issue_type="missing_rrna_type",
            severity="warning",
            details=f"no {kind} rRNA annotation found",
        )
        for kind in _STANDARD_RRNA_TYPES
        if kind not in observed_rrna
    ]
    findings.extend(
        AnnotationIssue(
            issue_type="missing_trna_amino_acid",
            severity="warning",
            details=f"no tRNA annotation found for {amino_acid}",
        )
        for amino_acid in sorted(_STANDARD_TRNA_AMINO_ACIDS - observed_trna)
    )
    return findings


def feature_quality_findings(
    features: Iterable[Feature],
    genome: Mapping[str, str] | None = None,
    circular_seqids: frozenset[str] = frozenset(),
) -> list[AnnotationIssue]:
    """Return coordinate, pseudogene, partial, adjacency, and RNA findings."""

    feature_list = list(features)
    cds_features = [feature for feature in feature_list if feature.type == "CDS"]
    findings: list[AnnotationIssue] = []
    for feature in feature_list:
        if _is_pseudogene(feature):
            findings.append(
                _finding_for_feature(
                    "pseudogene",
                    "info",
                    feature,
                    "feature explicitly marked pseudo or pseudogene",
                )
            )
        if _is_partial(feature):
            findings.append(
                _finding_for_feature(
                    "partial_feature",
                    "warning",
                    feature,
                    "feature has an incomplete boundary annotation",
                )
            )
    findings.extend(_overlap_findings(cds_features, genome, circular_seqids))
    findings.extend(_adjacent_duplicate_findings(cds_features))
    findings.extend(_structural_rna_findings(feature_list))
    return findings


def sequence_quality_findings(
    record: CdsSequence,
    genetic_code: int = 11,
    *,
    sequence_length: int | None = None,
    circular: bool = False,
) -> list[AnnotationIssue]:
    """Return conservative sequence-level review findings for one CDS.

    Explicit partial, pseudo, translational-exception, and ribosomal-slippage
    annotations suppress checks whose assumptions no longer hold.
    """

    if genetic_code not in _START_CODONS:
        raise ValueError(f"unsupported NCBI genetic code {genetic_code}; choose 4, 11, or 25")
    findings: list[AnnotationIssue] = []
    coding_sequence = record.coding_nucleotide.upper()
    partial = _is_partial(record.feature)
    five_prime_partial = _is_five_prime_partial(record.feature)
    three_prime_partial = _is_three_prime_partial(record.feature)
    pseudo = _is_pseudogene(record.feature)
    translation_exception = bool(record.feature.attributes.get("transl_except"))
    biological_exception = translation_exception or bool(
        record.feature.attributes.get("exception")
    )
    segments = record.segments or (record.feature,)
    near_left_boundary = (
        sequence_length is not None and min(segment.start for segment in segments) <= 3
    )
    near_right_boundary = (
        sequence_length is not None
        and max(segment.end for segment in segments) >= sequence_length - 2
    )
    five_prime_boundary = not circular and (
        near_left_boundary if record.feature.strand != "-" else near_right_boundary
    )
    three_prime_boundary = not circular and (
        near_right_boundary if record.feature.strand != "-" else near_left_boundary
    )
    if len(coding_sequence) % 3 and not (
        partial or pseudo or biological_exception or five_prime_boundary or three_prime_boundary
    ):
        findings.append(
            _finding_for_feature(
                "non_triplet_cds",
                "warning",
                record.feature,
                f"coding length {len(coding_sequence)} is not divisible by three",
            )
        )
    exception_indices = set(record.translation_exception_indices)
    unexplained_internal_stop = any(
        amino_acid == "*" and index not in exception_indices
        for index, amino_acid in enumerate(record.protein)
    )
    if unexplained_internal_stop and not pseudo:
        findings.append(
            _finding_for_feature(
                "internal_stop_codon",
                "warning",
                record.feature,
                "translated CDS contains an internal stop codon",
            )
        )
    if (
        len(coding_sequence) >= 3
        and coding_sequence[:3] not in _START_CODONS[genetic_code]
        and not (five_prime_partial or pseudo or biological_exception or five_prime_boundary)
    ):
        findings.append(
            _finding_for_feature(
                "unrecognized_start_codon",
                "warning",
                record.feature,
                f"complete CDS begins with {coding_sequence[:3]}",
            )
        )
    if (
        len(coding_sequence) >= 3
        and len(coding_sequence) % 3 == 0
        and coding_sequence[-3:] not in _STOP_CODONS[genetic_code]
        and not (three_prime_partial or pseudo or biological_exception or three_prime_boundary)
    ):
        findings.append(
            _finding_for_feature(
                "missing_stop_codon",
                "warning",
                record.feature,
                f"complete CDS ends with {coding_sequence[-3:]}, not a stop codon",
            )
        )
    ambiguous_bases = sum(base not in "ACGT" for base in coding_sequence)
    if ambiguous_bases:
        findings.append(
            _finding_for_feature(
                "ambiguous_cds_bases",
                "warning",
                record.feature,
                f"coding sequence contains {ambiguous_bases} ambiguous bases",
            )
        )
    return findings


def _covered_bases(
    features: Iterable[Feature],
    genome: Mapping[str, str],
    circular_seqids: frozenset[str],
) -> int:
    """Return the union length of CDS intervals across all sequence records."""

    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for feature in features:
        if feature.type != "CDS" or feature.seqid not in genome:
            continue
        sequence_length = len(genome[feature.seqid])
        if feature.seqid in circular_seqids:
            intervals[feature.seqid].extend(
                _circular_intervals(feature, sequence_length)
            )
        elif feature.end <= sequence_length:
            intervals[feature.seqid].append((feature.start, feature.end))

    covered = 0
    for sequence_intervals in intervals.values():
        merged_start = merged_end = None
        for start, end in sorted(sequence_intervals):
            if merged_start is None:
                merged_start, merged_end = start, end
            elif start <= merged_end + 1:
                merged_end = max(merged_end, end)
            else:
                covered += merged_end - merged_start + 1
                merged_start, merged_end = start, end
        if merged_start is not None:
            covered += merged_end - merged_start + 1
    return covered


def quality_summary(
    features: Iterable[Feature],
    genome: Mapping[str, str],
    circular_seqids: frozenset[str],
    findings: Iterable[AnnotationIssue],
) -> dict[str, object]:
    """Summarize genome composition, coding density, and finding counts."""

    feature_list = list(features)
    finding_list = list(findings)
    genome_length = sum(map(len, genome.values()))
    gc_bases = sum(sequence.count("G") + sequence.count("C") for sequence in genome.values())
    explicit_cds_ids: set[str] = set()
    anonymous_cds_count = 0
    for feature in feature_list:
        if feature.type != "CDS":
            continue
        explicit_id = feature.attributes.get("ID")
        if explicit_id:
            explicit_cds_ids.add(explicit_id)
        else:
            anonymous_cds_count += 1
    cds_count = len(explicit_cds_ids) + anonymous_cds_count
    covered_bases = _covered_bases(feature_list, genome, circular_seqids)
    issue_counts = Counter(finding.issue_type for finding in finding_list)
    severity_counts = Counter(finding.severity for finding in finding_list)
    return {
        "genome_gc_percent": 100 * gc_bases / genome_length if genome_length else 0,
        "coding_density_percent": 100 * covered_bases / genome_length if genome_length else 0,
        "cds_per_kb": 1000 * cds_count / genome_length if genome_length else 0,
        "finding_count": len(finding_list),
        "issue_counts": dict(sorted(issue_counts.items())),
        "severity_counts": dict(sorted(severity_counts.items())),
    }
