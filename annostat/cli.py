"""Command-line entry point and workflow orchestration for annostat."""

from __future__ import annotations

import sys
import tracemalloc
from collections import Counter
from pathlib import Path
from time import perf_counter
from typing import Callable

# Direct execution places ``annostat/`` rather than the repository root on
# sys.path. Add the parent before importing the package modules below.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse

from annostat import __version__
from annostat.analysis import COG_CATEGORY_NAMES, analyze_features
from annostat.output import (
    write_cds_fastas,
    write_codon_usage,
    write_count_table,
    write_overview,
    write_summary,
)
from annostat.parsers import parse_fasta, parse_gff
from annostat.plots import write_bar_chart, write_histogram
from annostat.report import render_html_report
from annostat.sequences import iter_cds_sequences


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""

    parser = argparse.ArgumentParser(
        prog="annostat",
        description="Analyze bacterial GFF3 annotations against a genome FASTA file.",
        epilog=(
            "examples:\n"
            "  %(prog)s -f genome.fna -g annotations.gff3\n"
            "  %(prog)s -f genome.fna -g annotations.gff3 -o results --table-format tsv\n\n"
            "The output includes an offline HTML report, analysis tables, CDS FASTA files, "
            "and three publication-ready SVG charts."
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
    parser.add_argument(
        "--profile", action="store_true",
        help="measure and report peak Python memory in addition to stage timings",
    )
    parser.add_argument("--version", action="version", version=f"annostat {__version__}")
    return parser


def run_analysis(
    fasta_path: Path,
    gff_path: Path,
    output_dir: Path,
    table_format: str,
    progress: Callable[[str], None] | None = None,
    profile: bool = False,
) -> dict[str, object]:
    """Run the complete annotation-analysis workflow and return its summary.

    The workflow parses both inputs, streams CDS extraction into FASTA output,
    calculates annotation statistics, writes tables and plots, and finally builds
    the offline HTML report. Optional callbacks receive the five stage messages.
    """

    notify = progress or (lambda message: None)
    stage_timings: dict[str, float] = {}
    if profile:
        tracemalloc.start()

    notify("Reading GFF3 annotations and FASTA sequences")
    stage_started = perf_counter()
    features = list(parse_gff(gff_path))
    genome = parse_fasta(fasta_path)
    circular_seqids = frozenset(
        feature.seqid
        for feature in features
        if feature.type == "region"
        and feature.attributes.get("Is_circular", "").lower() == "true"
    )
    stage_timings["input_parsing"] = perf_counter() - stage_started

    output_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = output_dir / "tables"
    sequences_dir = output_dir / "sequences"
    plots_dir = output_dir / "plots"
    for directory in (tables_dir, sequences_dir, plots_dir):
        directory.mkdir(exist_ok=True)

    notify("Streaming, translating, and writing CDS sequences")
    stage_started = perf_counter()
    codon_counts: Counter[str] = Counter()
    start_counts: Counter[str] = Counter()
    cds_lengths: list[int] = []

    def observed_records():
        """Yield CDS records while collecting lengths for the histogram."""

        # FASTA writing consumes this generator, fusing extraction and export.
        for record in iter_cds_sequences(
            features,
            genome,
            circular_seqids,
            codon_counts=codon_counts,
            start_counts=start_counts,
        ):
            cds_lengths.append(record.feature.length)
            yield record

    write_cds_fastas(sequences_dir, observed_records())
    stage_timings["cds_processing"] = perf_counter() - stage_started

    notify("Calculating statistics and writing analysis tables")
    stage_started = perf_counter()
    summary = analyze_features(features)
    summary["annostat_version"] = __version__
    summary["input_files"] = {"fasta": str(fasta_path), "gff3": str(gff_path)}
    summary["sequence_ids"] = sorted(genome)
    summary["circular_sequence_ids"] = sorted(circular_seqids)
    summary["genome_length"] = sum(map(len, genome.values()))
    summary["complete_codon_count"] = codon_counts.total()
    summary["start_codon_counts"] = dict(sorted(start_counts.items()))
    summary["top_codons"] = [
        {
            "codon": codon,
            "count": count,
            "percentage": 100 * count / codon_counts.total() if codon_counts else 0,
        }
        for codon, count in codon_counts.most_common(5)
    ]

    delimiter = "," if table_format == "csv" else "\t"
    write_overview(tables_dir / f"features.{table_format}", features, delimiter)
    write_codon_usage(tables_dir / "codon_usage.csv", codon_counts)
    write_count_table(tables_dir / "start_codons.csv", "start_codon", start_counts)
    write_count_table(
        tables_dir / "cog_categories.csv",
        "cog_category",
        summary["cog_category_counts"],
    )
    stage_timings["table_generation"] = perf_counter() - stage_started

    notify("Rendering scientific visualizations")
    stage_started = perf_counter()
    cog_plot_counts = {
        f"{category} - {COG_CATEGORY_NAMES.get(category, 'Unclassified')}": count
        for category, count in summary["cog_category_counts"].items()
    }
    write_bar_chart(
        plots_dir / "cog_categories.svg",
        "COG category distribution",
        cog_plot_counts,
        description="Functional category assignments; multi-category proteins contribute to each category",
        axis_label="COG assignments",
        percentage_total=sum(cog_plot_counts.values()),
    )
    write_histogram(
        plots_dir / "cds_lengths.svg",
        "CDS length distribution (nucleotides)",
        cds_lengths,
    )
    grouped_starts = {
        "ATG": start_counts.get("ATG", 0),
        "GTG": start_counts.get("GTG", 0),
        "TTG": start_counts.get("TTG", 0),
        "Other": sum(start_counts.values())
        - start_counts.get("ATG", 0)
        - start_counts.get("GTG", 0)
        - start_counts.get("TTG", 0),
    }
    # The chart stays readable while the CSV retains every observed start codon.
    write_bar_chart(
        plots_dir / "start_codons.svg",
        "Start codon usage",
        grouped_starts,
        description=f"Observed first codon across {len(cds_lengths):,} coding sequences",
        axis_label="Coding sequences",
        sort_by_value=False,
        percentage_total=len(cds_lengths),
    )
    stage_timings["plot_generation"] = perf_counter() - stage_started

    notify("Building the offline HTML report")
    plot_paths = [
        (plots_dir / "cog_categories.svg", "COG functional categories"),
        (plots_dir / "cds_lengths.svg", "CDS length distribution"),
        (plots_dir / "start_codons.svg", "Start codon usage"),
    ]
    summary["output_files"] = [
        "report.html",
        "summary.json",
        f"tables/features.{table_format}",
        "tables/codon_usage.csv",
        "tables/start_codons.csv",
        "tables/cog_categories.csv",
        "sequences/cds_nucleotide.fasta",
        "sequences/cds_protein.fasta",
        "plots/cog_categories.svg",
        "plots/cds_lengths.svg",
        "plots/start_codons.svg",
    ]
    summary["performance"] = {
        "stage_seconds": dict(stage_timings),
        "total_seconds": sum(stage_timings.values()),
        "peak_memory_bytes": tracemalloc.get_traced_memory()[1] if profile else None,
    }
    stage_started = perf_counter()
    report_html = render_html_report(summary, plot_paths)
    if profile:
        summary["performance"]["peak_memory_bytes"] = tracemalloc.get_traced_memory()[1]
        report_html = render_html_report(summary, plot_paths)
    (output_dir / "report.html").write_text(
        report_html, encoding="utf-8"
    )
    stage_timings["report_generation"] = perf_counter() - stage_started
    peak_memory = tracemalloc.get_traced_memory()[1] if profile else None
    if profile:
        tracemalloc.stop()
    summary["performance"] = {
        "stage_seconds": dict(stage_timings),
        "total_seconds": sum(stage_timings.values()),
        "peak_memory_bytes": peak_memory,
    }
    write_summary(output_dir / "summary.json", summary)
    return summary


def _percentage(part: int, whole: int) -> str:
    """Format a part-to-whole ratio as a percentage, including zero totals."""

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
    print(f"  Files written            {len(summary['output_files'])}")
    print(f"  Completed in             {elapsed:.2f} seconds")
    performance = summary["performance"]
    if performance.get("peak_memory_bytes"):
        print(f"  Peak Python memory       {performance['peak_memory_bytes'] / 1024 / 1024:.2f} MiB")
        print("\nPerformance profile")
        print("-" * 56)
        stage_labels = {
            "input_parsing": "Input parsing",
            "cds_processing": "CDS processing and FASTA",
            "table_generation": "Analysis tables",
            "plot_generation": "SVG plots",
            "report_generation": "HTML report",
        }
        for stage, seconds in performance["stage_seconds"].items():
            print(f"  {stage_labels.get(stage, stage):<30} {seconds:>12.4f} s")
        print("-" * 56)


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
        """Print one numbered workflow-stage notification."""

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
            args.profile,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if not args.quiet:
        _print_summary(summary, args.output, perf_counter() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
