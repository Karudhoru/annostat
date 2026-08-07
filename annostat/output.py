"""Writers for annostat tables, summaries, and FASTA output."""

from __future__ import annotations

import csv
import json
from collections import Counter
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import TextIO

from annostat.models import CdsSequence, Feature


def write_overview(path: Path, features: Iterable[Feature], delimiter: str) -> None:
    """Write the assignment's compact CSV or TSV feature overview.

    Missing ID, gene, and product attributes are deliberately written as empty
    fields instead of placeholder text.
    """

    fieldnames = [
        "ID", "Feature type", "Start position", "Stop position",
        "Gene name", "Product", "Strand",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for feature in features:
            writer.writerow(
                {
                    "ID": feature.attributes.get("ID", ""),
                    "Feature type": feature.type,
                    "Start position": feature.start,
                    "Stop position": feature.end,
                    "Gene name": feature.attributes.get("gene", ""),
                    "Product": feature.attributes.get("product", ""),
                    "Strand": feature.strand,
                }
            )


def _write_fasta_record(handle: TextIO, header: str, sequence: str) -> None:
    """Write one FASTA record with sequence lines wrapped at 60 characters."""

    wrapped_sequence = "\n".join(
        sequence[offset : offset + 60] for offset in range(0, len(sequence), 60)
    )
    handle.write(f">{header}\n{wrapped_sequence}\n")


def write_cds_fastas(output_dir: Path, records: Iterable[CdsSequence]) -> None:
    """Write paired nucleotide and amino-acid multi-FASTA files in one pass."""

    with (
        (output_dir / "cds_nucleotide.fasta").open("w", encoding="utf-8") as nucleotide_file,
        (output_dir / "cds_protein.fasta").open("w", encoding="utf-8") as protein_file,
    ):
        for record in records:
            description = record.feature.attributes.get("product", "")
            header = record.feature.id + (f" {description}" if description else "")
            _write_fasta_record(nucleotide_file, header, record.nucleotide)
            _write_fasta_record(protein_file, header, record.protein)


def write_count_table(path: Path, heading: str, counts: Mapping[str, int]) -> None:
    """Write sorted labels and integer counts as a two-column CSV table."""

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([heading, "count"])
        writer.writerows(sorted(counts.items()))


def write_codon_usage(path: Path, counts: Counter[str]) -> None:
    """Write codon counts and percentages over all complete CDS codons."""

    total = counts.total()
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["codon", "count", "percentage"])
        for codon, count in sorted(counts.items()):
            percentage = (100 * count / total) if total else 0
            writer.writerow([codon, count, f"{percentage:.6f}"])


def write_summary(path: Path, summary: Mapping[str, object]) -> None:
    """Serialize the complete analysis summary as indented JSON."""

    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
