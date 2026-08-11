"""Command-line entry point and workflow orchestration for annostat."""

from __future__ import annotations

import hashlib
import json
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
from annostat.cohort import build_cohort, write_cohort
from annostat.comparison import GenomeInput, run_comparison
from annostat.filtering import CdsFilter
from annostat.ncbi import fetch_genomes
from annostat.output import (
    write_annotation_findings,
    write_cds_fastas,
    write_codon_usage,
    write_count_table,
    write_overview,
    write_summary,
)
from annostat.parsers import parse_fasta, parse_gff
from annostat.plots import write_bar_chart, write_histogram
from annostat.qc import (
    feature_quality_findings,
    quality_summary,
    sequence_quality_findings,
)
from annostat.report import render_html_report
from annostat.sequences import (
    SUPPORTED_GENETIC_CODES,
    declared_genetic_codes,
    iter_cds_sequences,
    recognized_start_codons,
    resolve_genetic_code,
)
from annostat.validation import validate_annotation, write_validation


def _scientific_fingerprint(summary: dict[str, object]) -> str:
    """Hash scientific outputs while excluding paths and performance timings."""

    payload = {
        key: value
        for key, value in summary.items()
        if key not in {"input_files", "output_files", "performance", "scientific_fingerprint"}
    }
    validation = payload.get("validation")
    if isinstance(validation, dict):
        payload["validation"] = {
            key: value for key, value in validation.items() if key != "input_files"
        }
    return hashlib.sha256(
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode("utf-8")
    ).hexdigest()


def build_root_parser() -> argparse.ArgumentParser:
    """Build the command overview shown by ``annostat --help``."""

    return argparse.ArgumentParser(
        prog="annostat",
        description=(
            "Analyze bacterial GFF3 annotations, validate input integrity, and "
            "summarize or compare completed analyses."
        ),
        epilog=(
            "commands:\n"
            "  inspect     analyze one FASTA/GFF3 annotation pair\n"
            "  validate    check deterministic structural integrity\n"
            "  summarize   aggregate completed Annostat summaries\n"
            "  compare     compare two or more annotated genomes\n"
            "  fetch       download annotated NCBI assemblies\n\n"
            "Run 'annostat COMMAND --help' for command-specific options. "
            "Legacy 'annostat -f ... -g ...' inspection remains supported."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )


def build_parser(prog: str = "annostat") -> argparse.ArgumentParser:
    """Build the command-line argument parser."""

    parser = argparse.ArgumentParser(
        prog=prog,
        description="Analyze bacterial GFF3 annotations against a genome FASTA file.",
        epilog=(
            "examples:\n"
            "  %(prog)s -f genome.fna -g annotations.gff3\n"
            "  annostat validate -f genome.fna -g annotations.gff3\n"
            "  annostat summarize results/batch -o cohort\n"
            "  %(prog)s -f genome.fna -g annotations.gff3 -o results --table-format tsv\n"
            "  annostat compare --genome a a.fna a.gff3 --genome b b.fna b.gff3\n"
            "  annostat fetch GCF_000007145.1 -o ncbi_data\n\n"
            "The output includes an offline HTML report, analysis tables, CDS FASTA files, "
            "and up to three publication-ready SVG charts. Use the compare and fetch commands "
            "for multi-genome profiles and NCBI assembly downloads."
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
    parser.add_argument(
        "--genetic-code", type=int, choices=(4, 11, 25), default=None, metavar="TABLE",
        help=(
            "override the NCBI translation table; otherwise use GFF3 transl_table "
            "or fall back to 11"
        ),
    )
    filter_group = parser.add_argument_group(
        "filtered CDS export",
        "write additional matching CDS tables and FASTA files without changing the full analysis",
    )
    filter_group.add_argument(
        "--min-cds-length",
        type=_positive_int,
        metavar="BP",
        help="include CDS features at least BP nucleotides long",
    )
    filter_group.add_argument(
        "--max-cds-length",
        type=_positive_int,
        metavar="BP",
        help="include CDS features at most BP nucleotides long",
    )
    filter_group.add_argument(
        "--require-cog",
        action="store_true",
        help="include only CDS features with at least one COG category",
    )
    filter_group.add_argument(
        "--exclude-hypothetical",
        action="store_true",
        help="exclude CDS features annotated as hypothetical proteins",
    )
    parser.add_argument("--version", action="version", version=f"Annostat {__version__}")
    return parser


def build_validate_parser() -> argparse.ArgumentParser:
    """Build the parser for deterministic FASTA/GFF3 validation."""

    parser = argparse.ArgumentParser(
        prog="annostat validate",
        description=(
            "Validate structural and cross-file integrity without making "
            "taxon-dependent biological claims."
        ),
        epilog=(
            "example:\n"
            "  annostat validate -f genome.fna -g annotations.gff3 -o validation\n\n"
            "Exit status 0 means the selected failure threshold was not reached; "
            "status 1 means it was reached. JSON and TSV results are always written."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-f", "--fasta", required=True, type=Path, metavar="FILE",
        help="genome FASTA file",
    )
    parser.add_argument(
        "-g", "--gff", required=True, type=Path, metavar="FILE",
        help="matching GFF3 annotation file",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("annostat_validation"), metavar="DIR",
        help="result directory (default: annostat_validation)",
    )
    parser.add_argument(
        "--fail-on", choices=("error", "warning", "never"), default="error",
        help="return status 1 for this severity or above (default: error)",
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true",
        help="suppress the terminal validation summary",
    )
    return parser


def build_summarize_parser() -> argparse.ArgumentParser:
    """Build the parser for MultiQC-style aggregation of Annostat outputs."""

    parser = argparse.ArgumentParser(
        prog="annostat summarize",
        description="Aggregate one or more Annostat summary.json files into a cohort report.",
        epilog=(
            "example:\n"
            "  annostat summarize results/sample-* -o cohort\n\n"
            "Inputs may be individual summary.json files or directories scanned recursively."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "inputs", nargs="+", type=Path, metavar="PATH",
        help="summary.json file or directory to scan recursively",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("annostat_cohort"), metavar="DIR",
        help="result directory (default: annostat_cohort)",
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true",
        help="suppress the terminal cohort summary",
    )
    return parser


def build_compare_parser() -> argparse.ArgumentParser:
    """Build the parser for comparative annotation profiling."""

    parser = argparse.ArgumentParser(
        prog="annostat compare",
        description="Compare two or more bacterial FASTA/GFF3 annotation datasets.",
        epilog=(
            "examples:\n"
            "  annostat compare --genome genome-a a.fna a.gff3 --genome genome-b b.fna b.gff3\n"
            "  annostat compare --genome sample sample.fna sample.gff3 --reference GCF_000007145.1\n\n"
            "Local assemblies and external NCBI references may be combined. NCBI inputs require the official "
            "NCBI Datasets CLI."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--genome", action="append", nargs=3, metavar=("LABEL", "FASTA", "GFF3"),
        help="labelled local genome input; repeat for each dataset",
    )
    parser.add_argument(
        "--reference", "--accession", dest="references", action="append", metavar="GCF_OR_GCA",
        help="external NCBI assembly to include in the comparison; repeat as needed",
    )
    parser.add_argument(
        "--genetic-code",
        dest="genetic_codes",
        action="append",
        nargs=2,
        metavar=("LABEL", "TABLE"),
        help=(
            "override the translation table for one labelled genome or reference; "
            "repeat as needed"
        ),
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("annostat_comparison"), metavar="DIR",
        help="output directory (default: annostat_comparison)",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="suppress progress and summary output")
    return parser


def build_fetch_parser() -> argparse.ArgumentParser:
    """Build the parser for explicit NCBI assembly downloads."""

    parser = argparse.ArgumentParser(
        prog="annostat fetch",
        description="Download annotated assemblies through the official NCBI Datasets CLI.",
        epilog=(
            "example:\n"
            "  annostat fetch GCF_000007145.1 -o ncbi_data\n\n"
            "Accessions must be versioned GCF_ or GCA_ identifiers. Existing downloaded "
            "assemblies are reused."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "accessions", nargs="+", metavar="GCF_OR_GCA",
        help="one or more versioned NCBI assembly accessions",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("annostat_ncbi"), metavar="DIR",
        help="extraction directory (default: annostat_ncbi)",
    )
    return parser


def _positive_int(value: str) -> int:
    """Parse a strictly positive integer for a command-line option."""

    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def run_analysis(
    fasta_path: Path,
    gff_path: Path,
    output_dir: Path,
    table_format: str,
    progress: Callable[[str], None] | None = None,
    profile: bool = False,
    cds_filter: CdsFilter | None = None,
    genetic_code: int | None = None,
) -> dict[str, object]:
    """Run the complete annotation-analysis workflow and return its summary.

    The workflow parses both inputs, streams CDS extraction into FASTA output,
    calculates annotation statistics, writes tables and plots, and finally builds
    the offline HTML report. Optional callbacks receive the five stage messages.
    """

    notify = progress or (lambda message: None)
    active_filter = cds_filter or CdsFilter()
    stage_timings: dict[str, float] = {}
    if profile:
        tracemalloc.start()

    notify("Reading GFF3 annotations and FASTA sequences")
    stage_started = perf_counter()
    validation = validate_annotation(fasta_path, gff_path)
    if not validation["valid"]:
        errors = [
            finding for finding in validation["findings"]
            if finding["severity"] == "error"
        ]
        first = errors[0]
        raise ValueError(
            f"input validation failed with {len(errors)} error(s): "
            f"{first['rule_id']}: {first['message']}"
        )
    features = list(parse_gff(gff_path))
    genome = parse_fasta(fasta_path)
    declared_codes = declared_genetic_codes(features)
    selected_genetic_code = resolve_genetic_code(features, genetic_code)
    genetic_code_source = (
        "command_line" if genetic_code is not None
        else "gff3" if declared_codes
        else "default"
    )
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
    stale_table_format = "tsv" if table_format == "csv" else "csv"
    (tables_dir / f"features.{stale_table_format}").unlink(missing_ok=True)
    filtered_dir = output_dir / "filtered"
    if active_filter.active:
        filtered_dir.mkdir(exist_ok=True)
        (filtered_dir / f"features.{stale_table_format}").unlink(missing_ok=True)
    elif filtered_dir.is_dir():
        for generated_name in (
            "features.csv",
            "features.tsv",
            "cds_nucleotide.fasta",
            "cds_protein.fasta",
        ):
            (filtered_dir / generated_name).unlink(missing_ok=True)
        try:
            filtered_dir.rmdir()
        except OSError:
            # Preserve any files in filtered/ that Annostat did not generate.
            pass

    notify("Streaming, translating, and writing CDS sequences")
    stage_started = perf_counter()
    codon_counts: Counter[str] = Counter()
    start_counts: Counter[str] = Counter()
    cds_lengths: list[int] = []
    filtered_records = []
    quality_findings = feature_quality_findings(features, genome, circular_seqids)

    def observed_records():
        """Yield CDS records while collecting lengths for the histogram."""

        # FASTA writing consumes this generator, fusing extraction and export.
        for record in iter_cds_sequences(
            features,
            genome,
            circular_seqids,
            codon_counts=codon_counts,
            start_counts=start_counts,
            genetic_code=selected_genetic_code,
        ):
            cds_lengths.append(record.length)
            quality_findings.extend(
                sequence_quality_findings(
                    record,
                    selected_genetic_code,
                    sequence_length=len(genome[record.feature.seqid]),
                    circular=record.feature.seqid in circular_seqids,
                )
            )
            if active_filter.active and active_filter.matches(
                record.feature, length=record.length
            ):
                filtered_records.append(record)
            yield record

    write_cds_fastas(sequences_dir, observed_records())
    if active_filter.active:
        write_cds_fastas(filtered_dir, filtered_records)
    stage_timings["cds_processing"] = perf_counter() - stage_started

    notify("Calculating statistics and writing analysis tables")
    stage_started = perf_counter()
    summary = analyze_features(features)
    summary["cog_data_available"] = bool(summary["cog_category_counts"])
    summary["annostat_version"] = __version__
    summary["genetic_code"] = selected_genetic_code
    summary["genetic_code_source"] = genetic_code_source
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
    summary["quality_control"] = quality_summary(
        features,
        genome,
        circular_seqids,
        quality_findings,
    )
    summary["validation"] = validation
    summary["filtered_export"] = {
        **active_filter.as_dict(),
        "selected_cds_count": len(filtered_records),
    }

    delimiter = "," if table_format == "csv" else "\t"
    write_overview(tables_dir / f"features.{table_format}", features, delimiter)
    write_codon_usage(tables_dir / "codon_usage.csv", codon_counts)
    write_count_table(tables_dir / "start_codons.csv", "start_codon", start_counts)
    write_count_table(
        tables_dir / "cog_categories.csv",
        "cog_category",
        summary["cog_category_counts"],
    )
    write_annotation_findings(
        tables_dir / "annotation_issues.csv",
        quality_findings,
    )
    write_validation(output_dir / "validation", validation)
    if active_filter.active:
        write_overview(
            filtered_dir / f"features.{table_format}",
            (record.feature for record in filtered_records),
            delimiter,
        )
    stage_timings["table_generation"] = perf_counter() - stage_started

    notify("Rendering scientific visualizations")
    stage_started = perf_counter()
    cog_plot_counts = {
        f"{category} - {COG_CATEGORY_NAMES.get(category, 'Unclassified')}": count
        for category, count in summary["cog_category_counts"].items()
    }
    cog_plot = plots_dir / "cog_categories.svg"
    cog_plot.unlink(missing_ok=True)
    if summary["cog_data_available"]:
        write_bar_chart(
            cog_plot,
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
    selected_start_codons = recognized_start_codons(selected_genetic_code)
    common_start_codons = {"ATG", "GTG", "TTG"}
    other_recognized_starts = sum(
        start_counts.get(codon, 0)
        for codon in selected_start_codons - common_start_codons
    )
    recognized_start_total = sum(
        start_counts.get(codon, 0) for codon in selected_start_codons
    )
    grouped_starts = {
        "ATG": start_counts.get("ATG", 0),
        "GTG": start_counts.get("GTG", 0),
        "TTG": start_counts.get("TTG", 0),
        "Other recognized": other_recognized_starts,
        "Unrecognized": sum(start_counts.values()) - recognized_start_total,
    }
    # The chart stays readable while the CSV retains every observed start codon.
    write_bar_chart(
        plots_dir / "start_codons.svg",
        "Start codon usage",
        grouped_starts,
        description=(
            f"Observed first codon across {len(cds_lengths):,} coding sequences; "
            f"recognized initiators follow NCBI table {selected_genetic_code}"
        ),
        axis_label="Coding sequences",
        sort_by_value=False,
        percentage_total=len(cds_lengths),
    )
    stage_timings["plot_generation"] = perf_counter() - stage_started

    notify("Building the offline HTML report")
    plot_paths = [
        (plots_dir / "cds_lengths.svg", "CDS length distribution"),
        (plots_dir / "start_codons.svg", "Start codon usage"),
    ]
    if summary["cog_data_available"]:
        plot_paths.insert(0, (cog_plot, "COG functional categories"))
    summary["output_files"] = [
        "report.html",
        "summary.json",
        f"tables/features.{table_format}",
        "tables/codon_usage.csv",
        "tables/start_codons.csv",
        "tables/cog_categories.csv",
        "tables/annotation_issues.csv",
        "validation/validation.json",
        "validation/validation.tsv",
        "sequences/cds_nucleotide.fasta",
        "sequences/cds_protein.fasta",
        "plots/cds_lengths.svg",
        "plots/start_codons.svg",
    ]
    if summary["cog_data_available"]:
        summary["output_files"].insert(-2, "plots/cog_categories.svg")
    if active_filter.active:
        summary["output_files"].extend(
            (
                f"filtered/features.{table_format}",
                "filtered/cds_nucleotide.fasta",
                "filtered/cds_protein.fasta",
            )
        )
    summary["performance"] = {
        "stage_seconds": dict(stage_timings),
        "total_seconds": sum(stage_timings.values()),
        "peak_memory_bytes": tracemalloc.get_traced_memory()[1] if profile else None,
    }
    summary["scientific_fingerprint"] = _scientific_fingerprint(summary)
    stage_started = perf_counter()
    render_html_report(summary, plot_paths)
    stage_timings["report_generation"] = perf_counter() - stage_started
    peak_memory = tracemalloc.get_traced_memory()[1] if profile else None
    if profile:
        tracemalloc.stop()
    summary["performance"] = {
        "stage_seconds": dict(stage_timings),
        "total_seconds": sum(stage_timings.values()),
        "peak_memory_bytes": peak_memory,
    }
    summary["scientific_fingerprint"] = _scientific_fingerprint(summary)
    (output_dir / "report.html").write_text(
        render_html_report(summary, plot_paths), encoding="utf-8"
    )
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
    recognized_starts = sum(
        start_counts.get(codon, 0)
        for codon in recognized_start_codons(int(summary["genetic_code"]))
    )
    quality = summary["quality_control"]
    warning_count = quality["severity_counts"].get("warning", 0)
    information_count = quality["severity_counts"].get("info", 0)
    warning_label = "warning" if warning_count == 1 else "warnings"
    information_label = "finding" if information_count == 1 else "findings"
    metrics = (
        ("Genome size", f"{int(summary['genome_length']):,} bp"),
        ("Sequences", f"{len(summary['sequence_ids']):,} ({len(summary['circular_sequence_ids']):,} circular)"),
        ("Features", f"{int(summary['total_features']):,}"),
        ("CDS", f"{cds_count:,}"),
        ("RNA features", f"{rna_count:,}"),
        ("Genome GC", f"{quality['genome_gc_percent']:.2f}%"),
        ("Coding density", f"{quality['coding_density_percent']:.2f}%"),
        (
            "Hypothetical CDS",
            f"{int(summary['hypothetical_cds_count']):,} ({_percentage(int(summary['hypothetical_cds_count']), cds_count)})",
        ),
        (
            "COG-annotated CDS",
            (
                f"{int(summary['cds_with_cog_count']):,} ({_percentage(int(summary['cds_with_cog_count']), cds_count)})"
                if summary["cog_data_available"] else "not available"
            ),
        ),
        (
            "Recognized starts",
            f"{recognized_starts:,} ({_percentage(recognized_starts, cds_count)})",
        ),
        (
            "QC review",
            f"{warning_count:,} {warning_label}, "
            f"{information_count:,} informational {information_label}",
        ),
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

    command_arguments = list(sys.argv[1:] if arguments is None else arguments)
    if not command_arguments or command_arguments == ["--help"] or command_arguments == ["-h"]:
        build_root_parser().print_help()
        return 0
    if command_arguments == ["--version"]:
        print(f"Annostat {__version__}")
        return 0
    if command_arguments[:1] == ["validate"]:
        return _main_validate(command_arguments[1:])
    if command_arguments[:1] == ["inspect"]:
        return _main_analysis(command_arguments[1:], prog="annostat inspect")
    if command_arguments[:1] == ["summarize"]:
        return _main_summarize(command_arguments[1:])
    if command_arguments[:1] == ["compare"]:
        return _main_compare(command_arguments[1:])
    if command_arguments[:1] == ["fetch"]:
        return _main_fetch(command_arguments[1:])

    if command_arguments and not command_arguments[0].startswith("-"):
        build_root_parser().error(
            f"unknown command {command_arguments[0]!r}; choose inspect, validate, "
            "summarize, compare, or fetch"
        )

    return _main_analysis(command_arguments)


def _main_analysis(command_arguments: list[str], prog: str = "annostat") -> int:
    """Run the full biological inspection and reporting workflow."""

    parser = build_parser(prog)
    args = parser.parse_args(command_arguments)
    if (
        args.min_cds_length is not None
        and args.max_cds_length is not None
        and args.min_cds_length > args.max_cds_length
    ):
        parser.error("--min-cds-length cannot be greater than --max-cds-length")
    cds_filter = CdsFilter(
        min_length=args.min_cds_length,
        max_length=args.max_cds_length,
        require_cog=args.require_cog,
        exclude_hypothetical=args.exclude_hypothetical,
    )
    started = perf_counter()
    if not args.quiet:
        print(f"Annostat {__version__} | bacterial genome annotation analysis")
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
            profile=args.profile,
            cds_filter=cds_filter,
            genetic_code=args.genetic_code,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if not args.quiet:
        _print_summary(summary, args.output, perf_counter() - started)
    return 0


def _main_validate(arguments: list[str]) -> int:
    """Run deterministic structural validation and return a CI-friendly status."""

    parser = build_validate_parser()
    args = parser.parse_args(arguments)
    result = validate_annotation(args.fasta, args.gff)
    files = write_validation(args.output, result)
    counts = result["severity_counts"]
    if not args.quiet:
        status = "PASS" if result["valid"] else "FAIL"
        print(f"Annostat {__version__} | annotation validation {status}")
        print(f"  Errors                   {counts['error']:,}")
        print(f"  Warnings                 {counts['warning']:,}")
        print(f"  Findings                 {sum(counts.values()):,}")
        print(f"  Output                   {args.output.resolve()}")
        print(f"  Files written            {len(files)}")
    if args.fail_on == "never":
        return 0
    if counts["error"]:
        return 1
    if args.fail_on == "warning" and counts["warning"]:
        return 1
    return 0


def _main_summarize(arguments: list[str]) -> int:
    """Aggregate completed inspections into a deterministic cohort package."""

    parser = build_summarize_parser()
    args = parser.parse_args(arguments)
    try:
        cohort = build_cohort(args.inputs)
        files = write_cohort(args.output, cohort)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if not args.quiet:
        print(f"Annostat {__version__} | cohort annotation QC")
        print(f"  Samples                  {cohort['sample_count']:,}")
        print(f"  Output                   {args.output.resolve()}")
        print(f"  Files written            {len(files)}")
        print(f"  Report                   {(args.output / 'cohort.html').resolve()}")
    return 0


def _main_compare(arguments: list[str]) -> int:
    """Parse and run the comparative-analysis command."""

    parser = build_compare_parser()
    args = parser.parse_args(arguments)
    code_overrides: dict[str, int] = {}
    for label, raw_code in args.genetic_codes or []:
        if label in code_overrides:
            parser.error(f"--genetic-code was supplied more than once for {label!r}")
        try:
            code = int(raw_code)
        except ValueError:
            parser.error(f"genetic code for {label!r} must be 4, 11, or 25")
        if code not in SUPPORTED_GENETIC_CODES:
            parser.error(f"genetic code for {label!r} must be 4, 11, or 25")
        code_overrides[label] = code
    datasets = [
        GenomeInput(label, Path(fasta), Path(gff), genetic_code=code_overrides.get(label))
        for label, fasta, gff in (args.genome or [])
    ]
    try:
        requested_count = len(datasets) + len(args.references or [])
        if requested_count < 2:
            parser.error("provide at least two inputs using --genome and/or --reference")
        requested_labels = [dataset.label for dataset in datasets] + (args.references or [])
        if len(requested_labels) != len(set(requested_labels)):
            parser.error("comparison genome labels and references must be unique")
        unknown_overrides = sorted(set(code_overrides) - set(requested_labels))
        if unknown_overrides:
            parser.error(
                "--genetic-code label does not match an input: "
                + ", ".join(unknown_overrides)
            )
        if args.references:
            fetched = fetch_genomes(args.references, args.output / "external_inputs")
            datasets.extend(
                GenomeInput(
                    item.accession,
                    item.fasta,
                    item.gff,
                    item.metadata,
                    code_overrides.get(item.accession),
                )
                for item in fetched
            )
        if not args.quiet:
            print(f"Annostat {__version__} | comparative annotation analysis")
            for dataset in datasets:
                print(f"  {dataset.label}: {dataset.fasta} + {dataset.gff}")
            print()
        step = 0

        def report_progress(message: str) -> None:
            """Print one numbered comparison-stage notification."""

            nonlocal step
            step += 1
            print(f"[{step}/4] {message}")

        summary = run_comparison(
            datasets, args.output, None if args.quiet else report_progress
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if not args.quiet:
        print("\nComparison summary")
        print("-" * 60)
        print(f"  Datasets                 {len(summary['datasets']):>34,}")
        print(f"  Taxonomic scope          {summary['taxonomic_scope'].replace('_', ' '):>34}")
        print(f"  Pairwise comparisons     {len(summary['pairwise_comparisons']):>34,}")
        print(f"  Interpretation notes     {len(summary['warnings']):>34,}")
        print("-" * 60)
        print(f"  Report                   {(args.output / 'comparison.html').resolve()}")
    return 0


def _main_fetch(arguments: list[str]) -> int:
    """Parse and run an explicit NCBI Datasets download."""

    parser = build_fetch_parser()
    args = parser.parse_args(arguments)
    try:
        fetched = fetch_genomes(args.accessions, args.output)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Downloaded {len(fetched)} annotated NCBI genome(s) to {args.output.resolve()}")
    for item in fetched:
        organism = item.metadata.get("organism_name") or "organism not provided"
        print(f"  {item.accession} ({organism}): {item.fasta} + {item.gff}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
