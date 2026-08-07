"""Command-line entry point and workflow orchestration for annostat."""

from __future__ import annotations

import sys
from pathlib import Path
from time import perf_counter
from typing import Callable

# Direct execution places ``annostat/`` rather than the repository root on
# sys.path. Add the parent before importing the package modules below.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse

from annostat import __version__
from annostat.analysis import COG_CATEGORY_NAMES, analyze_codons, analyze_features
from annostat.output import (
    write_cds_fastas,
    write_codon_usage,
    write_count_table,
    write_overview,
    write_summary,
)
from annostat.parsers import parse_fasta, parse_gff
from annostat.plots import write_bar_chart, write_codon_heatmap, write_histogram
from annostat.sequences import extract_cds_sequences


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""

    parser = argparse.ArgumentParser(
        prog="annostat",
        description="Analyze bacterial GFF3 annotations against a genome FASTA file.",
        epilog=(
            "examples:\n"
            "  %(prog)s -f genome.fna -g annotations.gff3\n"
            "  %(prog)s -f genome.fna -g annotations.gff3 -o results --table-format tsv\n\n"
            "The output includes summary tables, CDS FASTA files, and publication-ready SVG charts."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-f", "--fasta", required=True, type=Path, metavar="FILE", help="genome FASTA file")
    parser.add_argument("-g", "--gff", required=True, type=Path, metavar="FILE", help="GFF3 annotation file")
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("annostat_output"), metavar="DIR",
        help="output directory (default: annostat_output)",
    )
    parser.add_argument(
        "--table-format", choices=("csv", "tsv"), default="csv",
        help="feature overview format (default: csv)",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="suppress progress and summary output")
    parser.add_argument("--version", action="version", version=f"annostat {__version__}")
    return parser


def run_analysis(
    fasta_path: Path,
    gff_path: Path,
    output_dir: Path,
    table_format: str,
    progress: Callable[[str], None] | None = None,
) -> dict[str, object]:
    """Run the complete annotation-analysis workflow."""

    notify = progress or (lambda message: None)
    notify("Reading GFF3 annotations and FASTA sequences")
    features = list(parse_gff(gff_path))
    genome = parse_fasta(fasta_path)
    circular_seqids = frozenset(
        feature.seqid
        for feature in features
        if feature.type == "region"
        and feature.attributes.get("Is_circular", "").lower() == "true"
    )
    notify("Extracting and translating CDS sequences")
    records = extract_cds_sequences(features, genome, circular_seqids)
    notify("Calculating feature, COG, and codon statistics")
    summary = analyze_features(features)
    codon_counts, start_counts = analyze_codons(records)
    summary["annostat_version"] = __version__
    summary["input_files"] = {"fasta": str(fasta_path), "gff3": str(gff_path)}
    summary["sequence_ids"] = sorted(genome)
    summary["circular_sequence_ids"] = sorted(circular_seqids)
    summary["genome_length"] = sum(map(len, genome.values()))
    summary["complete_codon_count"] = codon_counts.total()
    summary["start_codon_counts"] = dict(sorted(start_counts.items()))

    notify("Writing tables, summaries, and FASTA files")
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

    notify("Rendering scientific visualizations")
    cog_plot_counts = {
        f"{category} - {COG_CATEGORY_NAMES.get(category, 'Unclassified')}": count
        for category, count in summary["cog_category_counts"].items()
    }
    write_bar_chart(
        output_dir / "cog_categories.svg",
        "COG category distribution",
        cog_plot_counts,
        description="Functional category assignments; multi-category proteins contribute to each category",
        axis_label="COG assignments",
    )
    write_histogram(
        output_dir / "cds_lengths.svg",
        "CDS length distribution (nucleotides)",
        [record.feature.length for record in records],
    )
    write_bar_chart(
        output_dir / "start_codons.svg",
        "Start codon usage",
        start_counts,
        description=f"Observed first codon across {len(records):,} coding sequences",
        axis_label="Coding sequences",
    )
    write_codon_heatmap(output_dir / "codon_usage.svg", codon_counts)
    return summary


def _percentage(part: int, whole: int) -> str:
    return f"{100 * part / whole:.2f}%" if whole else "0.00%"


def _print_summary(summary: dict[str, object], output_dir: Path, elapsed: float) -> None:
    """Print a compact, human-readable analysis report."""

    cds_count = int(summary["cds_count"])
    rna_count = sum(summary["rna_counts"].values())
    start_counts = summary["start_codon_counts"]
    standard_starts = sum(start_counts.get(codon, 0) for codon in ("ATG", "GTG", "TTG"))
    metrics = (
        ("Genome size", f"{int(summary['genome_length']):,} bp"),
        ("Sequences", f"{len(summary['sequence_ids']):,} ({len(summary['circular_sequence_ids']):,} circular)"),
        ("Features", f"{int(summary['total_features']):,}"),
        ("CDS", f"{cds_count:,}"),
        ("RNA features", f"{rna_count:,}"),
        (
            "Hypothetical CDS",
            f"{int(summary['hypothetical_cds_count']):,} ({_percentage(int(summary['hypothetical_cds_count']), cds_count)})",
        ),
        (
            "COG-annotated CDS",
            f"{int(summary['cds_with_cog_count']):,} ({_percentage(int(summary['cds_with_cog_count']), cds_count)})",
        ),
        ("ATG/GTG/TTG starts", f"{standard_starts:,} ({_percentage(standard_starts, cds_count)})"),
    )
    print("\nAnalysis summary")
    print("-" * 56)
    for label, value in metrics:
        print(f"  {label:<24} {value:>28}")
    print("-" * 56)
    print(f"  Output                   {output_dir.resolve()}")
    print(f"  Files written            {len(list(output_dir.iterdir()))}")
    print(f"  Completed in             {elapsed:.2f} seconds")


def main(arguments: list[str] | None = None) -> int:
    """Parse arguments, run the analysis, and report its output location."""

    parser = build_parser()
    args = parser.parse_args(arguments)
    started = perf_counter()
    if not args.quiet:
        print(f"annostat {__version__} | bacterial genome annotation analysis")
        print(f"FASTA: {args.fasta.resolve()}")
        print(f"GFF3:  {args.gff.resolve()}\n")
    step = 0

    def report_progress(message: str) -> None:
        nonlocal step
        step += 1
        print(f"[{step}/5] {message}")

    try:
        summary = run_analysis(
            args.fasta,
            args.gff,
            args.output,
            args.table_format,
            None if args.quiet else report_progress,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if not args.quiet:
        _print_summary(summary, args.output, perf_counter() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
