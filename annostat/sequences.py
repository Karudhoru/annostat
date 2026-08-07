"""CDS extraction and translation utilities."""

from __future__ import annotations

from collections import Counter
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
_BACTERIAL_START_CODONS = {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}


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
) -> str:
    """Translate complete codons using NCBI bacterial genetic code 11.

    Unknown or ambiguous codons become ``X``. Valid alternative bacterial start
    codons translate to methionine in the first position, and a terminal stop is
    omitted from the returned protein sequence. If supplied, ``codon_counts`` is
    updated during the same traversal used for translation.
    """

    sequence = sequence.upper()
    codons = [
        sequence[offset : offset + 3] for offset in range(0, len(sequence) - 2, 3)
    ]
    if codon_counts is not None:
        codon_counts.update(codons)
    amino_acids = [_CODON_TABLE.get(codon, "X") for codon in codons]
    if codons and codons[0] in _BACTERIAL_START_CODONS:
        # Translation table 11 treats alternative starts as initiator methionine.
        amino_acids[0] = "M"
    if amino_acids and amino_acids[-1] == "*":
        amino_acids.pop()
    return "".join(amino_acids)


def extract_cds_sequences(
    features: Iterable[Feature],
    genome: Mapping[str, str],
    circular_seqids: frozenset[str] = frozenset(),
) -> list[CdsSequence]:
    """Return extracted and translated records for every CDS feature."""

    return list(iter_cds_sequences(features, genome, circular_seqids))


def iter_cds_sequences(
    features: Iterable[Feature],
    genome: Mapping[str, str],
    circular_seqids: frozenset[str] = frozenset(),
    *,
    codon_counts: Counter[str] | None = None,
    start_counts: MutableMapping[str, int] | None = None,
) -> Iterable[CdsSequence]:
    """Yield CDS records while optionally accumulating codon statistics.

    The iterator keeps sequence export memory-efficient: callers can consume and
    write each CDS immediately instead of retaining every translated sequence.
    GFF3 phase is applied after strand orientation and before translation.
    """

    for feature in features:
        if feature.type != "CDS":
            continue
        nucleotide = extract_feature_sequence_with_topology(feature, genome, circular_seqids)
        # Phase is the number of leading bases skipped to reach the next codon.
        phase = feature.phase or 0
        coding_nucleotide = nucleotide[phase:]
        if start_counts is not None and len(coding_nucleotide) >= 3:
            start = coding_nucleotide[:3].upper()
            start_counts[start] = start_counts.get(start, 0) + 1
        yield CdsSequence(
            feature=feature,
            nucleotide=nucleotide,
            coding_nucleotide=coding_nucleotide,
            protein=translate_dna(coding_nucleotide, codon_counts=codon_counts),
        )
