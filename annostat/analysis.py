"""Feature, codon, and COG analyses."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable

from annostat.models import CdsSequence, Feature


COG_CATEGORY_NAMES = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control and division",
    "E": "Amino acid metabolism and transport",
    "F": "Nucleotide metabolism and transport",
    "G": "Carbohydrate metabolism and transport",
    "H": "Coenzyme metabolism",
    "I": "Lipid metabolism",
    "J": "Translation and ribosome biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall and membrane biogenesis",
    "N": "Cell motility",
    "O": "Protein turnover and chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolite biosynthesis",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction",
    "U": "Intracellular trafficking and secretion",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages and transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}


def cog_categories(feature: Feature) -> tuple[str, ...]:
    """Return individual COG category letters from a feature's ``Dbxref``.

    A multi-category value such as ``COG:OE`` produces ``("O", "E")`` so each
    functional assignment contributes independently to the category totals.
    """

    categories: list[str] = []
    for cross_reference in feature.attributes.get("Dbxref", "").split(","):
        namespace, separator, value = cross_reference.partition(":")
        if separator and namespace == "COG" and value.isalpha() and value.isupper():
            # COG category strings contain one uppercase letter per assignment.
            categories.extend(value)
    return tuple(categories)


def analyze_features(features: Iterable[Feature]) -> dict[str, object]:
    """Calculate feature, RNA, annotation, and COG metrics in one pass.

    Counts that apply specifically to coding sequences are updated only for
    features whose GFF3 type is ``CDS``. The returned mapping is JSON-safe and
    can therefore be written directly into the analysis summary.
    """

    total_features = 0
    cds_count = 0
    hypothetical_cds_count = 0
    cds_with_gene_count = 0
    cds_with_cog_count = 0
    type_counts: Counter[str] = Counter()
    rna_counts: Counter[str] = Counter()
    cog_counts: Counter[str] = Counter()
    seen_features: set[tuple[str, str]] = set()
    gff_row_count = 0
    for row_index, feature in enumerate(features):
        gff_row_count += 1
        explicit_id = feature.attributes.get("ID")
        identity = (
            (feature.type, explicit_id)
            if explicit_id
            else ("__anonymous_row__", str(row_index))
        )
        if identity in seen_features:
            # Repeated IDs are rows of one multipart GFF3 feature.
            continue
        seen_features.add(identity)
        total_features += 1
        type_counts[feature.type] += 1
        if feature.type.lower().endswith("rna"):
            rna_counts[feature.type] += 1
        if feature.type != "CDS":
            continue
        cds_count += 1
        hypothetical_cds_count += (
            "hypothetical protein" in feature.attributes.get("product", "").lower()
        )
        cds_with_gene_count += "gene" in feature.attributes
        categories = cog_categories(feature)
        cds_with_cog_count += bool(categories)
        cog_counts.update(categories)

    return {
        "total_features": total_features,
        "gff_row_count": gff_row_count,
        "cds_count": cds_count,
        "rna_counts": dict(sorted(rna_counts.items())),
        "feature_type_counts": dict(sorted(type_counts.items())),
        "hypothetical_cds_count": hypothetical_cds_count,
        "cds_with_gene_count": cds_with_gene_count,
        "cds_with_cog_count": cds_with_cog_count,
        "cog_category_counts": dict(sorted(cog_counts.items())),
    }


def analyze_codons(records: Iterable[CdsSequence]) -> tuple[Counter[str], Counter[str]]:
    """Count pooled complete codons and observed start codons across CDS records.

    Incomplete trailing bases are excluded because they do not form a codon.
    The result is a pair containing all-codon counts followed by start counts.
    """

    codons: Counter[str] = Counter()
    starts: Counter[str] = Counter()
    for record in records:
        sequence = record.coding_nucleotide.upper()
        if len(sequence) >= 3:
            starts[sequence[:3]] += 1
        codons.update(sequence[offset : offset + 3] for offset in range(0, len(sequence) - 2, 3))
    return codons, starts
