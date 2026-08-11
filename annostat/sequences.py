"""CDS extraction and translation utilities."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, MutableMapping

from annostat.models import CdsSequence, Feature


_COMPLEMENT = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_GENETIC_CODE_OVERRIDES = {
    4: {"TGA": "W"},
    11: {},
    25: {"TGA": "G"},
}
SUPPORTED_GENETIC_CODES = frozenset(_GENETIC_CODE_OVERRIDES)
_START_CODONS = {
    4: {"TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"},
    11: {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"},
    25: {"TTG", "ATG", "GTG"},
}
_TRANSL_EXCEPT_PATTERN = re.compile(
    r"pos:(?:complement\()?([0-9]+)\.\.([0-9]+)\)?,aa:([A-Za-z*]+)",
    re.IGNORECASE,
)
_EXCEPTION_AMINO_ACIDS = {
    "sec": "U",
    "pyl": "O",
    "term": "*",
    "other": "X",
    "*": "*",
}


def declared_genetic_codes(features: Iterable[Feature]) -> frozenset[int]:
    """Return supported NCBI translation tables declared by CDS or region rows.

    Malformed, unsupported, or conflicting declarations are rejected. The
    validation workflow reports the same conditions as structured findings;
    this guard prevents inspection from silently using the wrong table.
    """

    codes: set[int] = set()
    for feature in features:
        if feature.type not in {"CDS", "region"} or "transl_table" not in feature.attributes:
            continue
        raw_code = feature.attributes["transl_table"].strip()
        try:
            code = int(raw_code)
        except ValueError as error:
            raise ValueError(
                f"feature {feature.id!r} has invalid transl_table {raw_code!r}"
            ) from error
        if code not in SUPPORTED_GENETIC_CODES:
            raise ValueError(
                f"feature {feature.id!r} declares unsupported NCBI genetic code {code}; "
                "choose 4, 11, or 25"
            )
        codes.add(code)
    if len(codes) > 1:
        choices = ", ".join(map(str, sorted(codes)))
        raise ValueError(f"annotation features declare conflicting translation tables: {choices}")
    return frozenset(codes)


def resolve_genetic_code(
    features: Iterable[Feature],
    requested: int | None = None,
) -> int:
    """Resolve a safe translation table from GFF3 metadata and a CLI override.

    A declared ``transl_table`` takes precedence over the bacterial default of
    table 11. An explicit override is accepted only when it agrees with a table
    declared by the annotation.
    """

    if requested is not None and requested not in SUPPORTED_GENETIC_CODES:
        raise ValueError(
            f"unsupported NCBI genetic code {requested}; choose 4, 11, or 25"
        )
    declared = declared_genetic_codes(features)
    if declared:
        annotation_code = next(iter(declared))
        if requested is not None and requested != annotation_code:
            raise ValueError(
                f"requested genetic code {requested} conflicts with GFF3 transl_table "
                f"{annotation_code}"
            )
        return annotation_code
    return requested if requested is not None else 11


def recognized_start_codons(genetic_code: int) -> frozenset[str]:
    """Return the initiator codons defined for a supported translation table."""

    try:
        return frozenset(_START_CODONS[genetic_code])
    except KeyError as error:
        raise ValueError(
            f"unsupported NCBI genetic code {genetic_code}; choose 4, 11, or 25"
        ) from error


def reverse_complement(sequence: str) -> str:
    """Return an uppercase reverse complement, including IUPAC ambiguity codes."""

    return sequence.translate(_COMPLEMENT)[::-1].upper()


def extract_feature_sequence(feature: Feature, genome: Mapping[str, str]) -> str:
    """Extract a linear feature using GFF3's 1-based inclusive coordinates."""

    return extract_feature_sequence_with_topology(feature, genome, frozenset())


def extract_feature_sequence_with_topology(
    feature: Feature,
    genome: Mapping[str, str],
    circular_seqids: frozenset[str],
) -> str:
    """Extract and orient a feature, including circular-origin wrapping.

    Linear features extending beyond their reference sequence are rejected.
    Features on the minus strand are reverse-complemented after all referenced
    nucleotide segments have been joined.
    """

    try:
        reference = genome[feature.seqid]
    except KeyError as error:
        raise ValueError(
            f"GFF3 sequence {feature.seqid!r} is missing from the FASTA file"
        ) from error
    if feature.end > len(reference) and feature.seqid not in circular_seqids:
        raise ValueError(
            f"feature {feature.id!r} ends at {feature.end}, beyond "
            f"{feature.seqid!r} length {len(reference)}"
        )
    if feature.end <= len(reference):
        # Python's zero-based, end-exclusive slice needs only the start adjusted.
        sequence = reference[feature.start - 1 : feature.end]
    else:
        # Circular features may require more than one segment around the origin.
        remaining = feature.length
        position = (feature.start - 1) % len(reference)
        parts: list[str] = []
        while remaining:
            chunk_length = min(remaining, len(reference) - position)
            parts.append(reference[position : position + chunk_length])
            remaining -= chunk_length
            position = 0
        sequence = "".join(parts)
    return reverse_complement(sequence) if feature.strand == "-" else sequence.upper()


def translate_dna(
    sequence: str,
    *,
    codon_counts: Counter[str] | None = None,
    genetic_code: int = 11,
    complete_start: bool = True,
) -> str:
    """Translate complete codons using a supported prokaryotic NCBI code.

    Unknown or ambiguous codons become ``X``. Valid alternative bacterial start
    codons translate to methionine in the first position, and a terminal stop is
    omitted from the returned protein sequence. If supplied, ``codon_counts`` is
    updated during the same traversal used for translation.
    """

    if genetic_code not in _GENETIC_CODE_OVERRIDES:
        raise ValueError(f"unsupported NCBI genetic code {genetic_code}; choose 4, 11, or 25")
    sequence = sequence.upper()
    codon_table = {**_CODON_TABLE, **_GENETIC_CODE_OVERRIDES[genetic_code]}
    codons = [
        sequence[offset : offset + 3] for offset in range(0, len(sequence) - 2, 3)
    ]
    if codon_counts is not None:
        codon_counts.update(codons)
    amino_acids = [codon_table.get(codon, "X") for codon in codons]
    if complete_start and codons and codons[0] in _START_CODONS[genetic_code]:
        # NCBI translation tables define initiators separately from elongation.
        amino_acids[0] = "M"
    if amino_acids and amino_acids[-1] == "*":
        amino_acids.pop()
    return "".join(amino_acids)


def _apply_translation_exceptions(
    protein: str,
    parts: list[Feature],
    phase: int,
    sequence_length: int,
) -> tuple[str, tuple[int, ...]]:
    """Apply position-specific NCBI ``transl_except`` annotations.

    The coordinate map is built in biological order, so the same logic handles
    plus, minus, multipart, and circular CDS features.
    """

    exception_text = ",".join(
        part.attributes.get("transl_except", "") for part in parts
        if part.attributes.get("transl_except")
    )
    if not exception_text:
        return protein, ()
    coordinates: list[int] = []
    for part in parts:
        raw_positions = (
            range(part.start, part.end + 1)
            if part.strand != "-"
            else range(part.end, part.start - 1, -1)
        )
        coordinates.extend((position - 1) % sequence_length + 1 for position in raw_positions)
    coding_coordinates = coordinates[phase:]
    amino_acids = list(protein)
    applied: list[int] = []
    for match in _TRANSL_EXCEPT_PATTERN.finditer(exception_text):
        start, end = int(match.group(1)), int(match.group(2))
        replacement = _EXCEPTION_AMINO_ACIDS.get(match.group(3).lower())
        if replacement is None:
            continue
        target = {
            (position - 1) % sequence_length + 1
            for position in range(min(start, end), max(start, end) + 1)
        }
        for codon_index in range(len(coding_coordinates) // 3):
            codon_positions = set(coding_coordinates[codon_index * 3 : codon_index * 3 + 3])
            if codon_positions != target:
                continue
            if codon_index < len(amino_acids):
                amino_acids[codon_index] = replacement
                applied.append(codon_index)
                break
            # ``translate_dna`` removes a terminal stop. Restore it when an
            # explicit terminal exception recodes that codon as an amino acid.
            if codon_index == len(amino_acids) and replacement != "*":
                amino_acids.append(replacement)
                applied.append(codon_index)
                break
    return "".join(amino_acids), tuple(sorted(set(applied)))


def extract_cds_sequences(
    features: Iterable[Feature],
    genome: Mapping[str, str],
    circular_seqids: frozenset[str] = frozenset(),
    *,
    genetic_code: int = 11,
) -> list[CdsSequence]:
    """Return extracted and translated records for every CDS feature."""

    return list(
        iter_cds_sequences(features, genome, circular_seqids, genetic_code=genetic_code)
    )


def iter_cds_sequences(
    features: Iterable[Feature],
    genome: Mapping[str, str],
    circular_seqids: frozenset[str] = frozenset(),
    *,
    codon_counts: Counter[str] | None = None,
    start_counts: MutableMapping[str, int] | None = None,
    genetic_code: int = 11,
) -> Iterable[CdsSequence]:
    """Yield CDS records while optionally accumulating codon statistics.

    The iterator keeps sequence export memory-efficient: callers can consume and
    write each CDS immediately instead of retaining every translated sequence.
    GFF3 phase is applied after strand orientation and before translation.
    """

    grouped: dict[str, list[Feature]] = defaultdict(list)
    order: list[str] = []
    for row_number, feature in enumerate(features):
        if feature.type != "CDS":
            continue
        # Only explicit IDs join rows. Coordinate-based fallback IDs represent
        # independent anonymous features and must never merge accidentally.
        group_id = feature.attributes.get("ID") or f"__anonymous_{row_number}"
        if group_id not in grouped:
            order.append(group_id)
        grouped[group_id].append(feature)

    for group_id in order:
        parts = grouped[group_id]
        signatures = {(part.seqid, part.strand) for part in parts}
        if len(signatures) != 1:
            raise ValueError(
                f"multipart CDS {group_id!r} has inconsistent sequence IDs or strands"
            )
        ordered = sorted(
            parts,
            key=lambda item: item.start,
            reverse=parts[0].strand == "-",
        )
        feature = ordered[0]
        nucleotide = "".join(
            extract_feature_sequence_with_topology(part, genome, circular_seqids)
            for part in ordered
        )
        # Only the biologically first segment's phase trims the joined CDS.
        # Later phases describe frame continuity at segment boundaries.
        phase = feature.phase or 0
        coding_nucleotide = nucleotide[phase:]
        if start_counts is not None and len(coding_nucleotide) >= 3:
            start = coding_nucleotide[:3].upper()
            start_counts[start] = start_counts.get(start, 0) + 1
        attributes = feature.attributes
        boundary_marker = "end_range" if feature.strand == "-" else "start_range"
        generic_partial = (
            attributes.get("partial", "").strip().lower() in {"1", "true", "yes"}
            and "start_range" not in attributes
            and "end_range" not in attributes
        )
        protein = translate_dna(
            coding_nucleotide,
            codon_counts=codon_counts,
            genetic_code=genetic_code,
            complete_start=boundary_marker not in attributes and not generic_partial,
        )
        protein, exception_indices = _apply_translation_exceptions(
            protein,
            ordered,
            phase,
            len(genome[feature.seqid]),
        )
        yield CdsSequence(
            feature=feature,
            nucleotide=nucleotide,
            coding_nucleotide=coding_nucleotide,
            protein=protein,
            segments=tuple(ordered),
            translation_exception_indices=exception_indices,
        )
