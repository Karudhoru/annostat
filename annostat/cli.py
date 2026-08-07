"""Command-line entry point and workflow orchestration for annostat."""

from __future__ import annotations

import sys
from pathlib import Path

# Direct execution places ``annostat/`` rather than the repository root on
# sys.path. Add the parent before importing the package modules below.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse

from annostat.analysis import analyze_codons, analyze_features
from annostat.output import (
    write_cds_fastas,
    write_codon_usage,
    write_count_table,
    write_overview,
    write_summary,
)
from annostat.parsers import parse_fasta, parse_gff
from annostat.plots import write_bar_chart, write_histogram
from annostat.sequences import extract_cds_sequences


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""

    parser = argparse.ArgumentParser(
        description="Analyze bacterial GFF3 annotations against a genome FASTA file."
    )
    parser.add_argument("-f", "--fasta", required=True, type=Path, help="genome FASTA file")
    parser.add_argument("-g", "--gff", required=True, type=Path, help="GFF3 annotation file")
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("annostat_output"), help="output directory"
    )
    parser.add_argument(
        "--table-format", choices=("csv", "tsv"), default="csv", help="overview table format"
    )
    return parser


def run_analysis(
    fasta_path: Path, gff_path: Path, output_dir: Path, table_format: str
) -> dict[str, object]:
    """Run the complete annotation-analysis workflow."""

    features = list(parse_gff(gff_path))
    genome = parse_fasta(fasta_path)
    circular_seqids = frozenset(
        feature.seqid
        for feature in features
        if feature.type == "region"
        and feature.attributes.get("Is_circular", "").lower() == "true"
    )
    records = extract_cds_sequences(features, genome, circular_seqids)
    summary = analyze_features(features)
    codon_counts, start_counts = analyze_codons(records)
    summary["sequence_ids"] = sorted(genome)
    summary["circular_sequence_ids"] = sorted(circular_seqids)
    summary["genome_length"] = sum(map(len, genome.values()))
    summary["complete_codon_count"] = codon_counts.total()
    summary["start_codon_counts"] = dict(sorted(start_counts.items()))

    output_dir.mkdir(parents=True, exist_ok=True)
    delimiter = "," if table_format == "csv" else "\t"
    write_overview(output_dir / f"features.{table_format}", features, delimiter)
    write_cds_fastas(output_dir, records)
    write_codon_usage(output_dir / "codon_usage.csv", codon_counts)
    write_count_table(output_dir / "start_codons.csv", "start_codon", start_counts)
    write_count_table(
        output_dir / "cog_categories.csv",
        "cog_category",
        summary["cog_category_counts"],
    )
    write_summary(output_dir / "summary.json", summary)
    write_bar_chart(
        output_dir / "cog_categories.svg",
        "COG category distribution",
        summary["cog_category_counts"],
    )
    write_histogram(
        output_dir / "cds_lengths.svg",
        "CDS length distribution (nucleotides)",
        [record.feature.length for record in records],
    )
    return summary


def main(arguments: list[str] | None = None) -> int:
    """Parse arguments, run the analysis, and report its output location."""

    parser = build_parser()
    args = parser.parse_args(arguments)
    try:
        summary = run_analysis(args.fasta, args.gff, args.output, args.table_format)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(
        f"Analyzed {summary['cds_count']} CDS and {summary['total_features']} total features. "
        f"Results: {args.output.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
