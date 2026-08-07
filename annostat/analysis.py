"""Feature, codon, and COG analyses."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable

from annostat.models import CdsSequence, Feature


def cog_categories(feature: Feature) -> tuple[str, ...]:
    """Extract individual COG category letters from a Dbxref attribute."""

    categories: list[str] = []
    for cross_reference in feature.attributes.get("Dbxref", "").split(","):
        namespace, separator, value = cross_reference.partition(":")
        if separator and namespace == "COG" and value.isalpha() and value.isupper():
            categories.extend(value)
    return tuple(categories)


def analyze_features(features: Iterable[Feature]) -> dict[str, object]:
    """Calculate required GFF3 feature and annotation metrics."""

    feature_list = list(features)
    type_counts = Counter(feature.type for feature in feature_list)
    rna_counts = Counter(
        feature.type for feature in feature_list if feature.type.lower().endswith("rna")
    )
    cds_features = [feature for feature in feature_list if feature.type == "CDS"]
    cog_counts: Counter[str] = Counter()
    for feature in cds_features:
        cog_counts.update(cog_categories(feature))

    return {
        "total_features": len(feature_list),
        "cds_count": len(cds_features),
        "rna_counts": dict(sorted(rna_counts.items())),
        "feature_type_counts": dict(sorted(type_counts.items())),
        "hypothetical_cds_count": sum(
            "hypothetical protein" in feature.attributes.get("product", "").lower()
            for feature in cds_features
        ),
        "cds_with_gene_count": sum("gene" in feature.attributes for feature in cds_features),
        "cds_with_cog_count": sum(bool(cog_categories(feature)) for feature in cds_features),
        "cog_category_counts": dict(sorted(cog_counts.items())),
    }


def analyze_codons(records: Iterable[CdsSequence]) -> tuple[Counter[str], Counter[str]]:
    """Count all complete codons and the first codon of each CDS."""

    codons: Counter[str] = Counter()
    starts: Counter[str] = Counter()
    for record in records:
        sequence = record.coding_nucleotide.upper()
        if len(sequence) >= 3:
            starts[sequence[:3]] += 1
        codons.update(sequence[offset : offset + 3] for offset in range(0, len(sequence) - 2, 3))
    return codons, starts
