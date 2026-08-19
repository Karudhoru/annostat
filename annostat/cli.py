"""Command-line entry point for annostat."""

from __future__ import annotations

import sys
from pathlib import Path
from time import perf_counter

# Direct execution places ``annostat/`` rather than the repository root on
# sys.path. Add the parent before importing the package modules below.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse

from annostat import __version__
from annostat.cohort import build_cohort, write_cohort
from annostat.comparison import GenomeInput, run_comparison
from annostat.filtering import CdsFilter
from annostat.ncbi import fetch_genomes
from annostat.sequences import SUPPORTED_GENETIC_CODES, recognized_start_codons
from annostat.validation import validate_annotation, write_validation
from annostat.workflow import run_analysis


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
            "  analyze     analyze one FASTA/GFF3 annotation pair\n"
            "  validate    check deterministic structural integrity\n"
            "  summarize   aggregate completed Annostat summaries\n"
            "  compare     compare two or more annotated genomes\n"
            "  fetch       download annotated NCBI assemblies\n\n"
            "example:\n"
            "  annostat analyze -f genome.fna -g annotations.gff3 -o results\n\n"
            "Run 'annostat COMMAND --help' for command-specific options. "
            "'inspect' is a supported alias; the option-only form is deprecated."
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
            "  %(prog)s -f genome.fna -g annotations.gff3 -o results --table-format tsv\n"
            "  %(prog)s -f genome.fna -g annotations.gff3 --profile\n\n"
            "Writes an offline HTML report, JSON summary, analysis tables, CDS FASTA files, "
            "validation results, and publication-ready SVG charts. Exit status 0 indicates "
            "a completed analysis; invalid inputs or options return status 2."
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
            "status 1 means it was reached. JSON and TSV results are written only "
            "when validation finds an issue."
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
    if command_arguments[:1] == ["analyze"]:
        return _main_analysis(command_arguments[1:], prog="annostat analyze")
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
            f"unknown command {command_arguments[0]!r}; choose analyze, inspect, "
            "validate, summarize, compare, or fetch"
        )

    print(
        "annostat: warning: the option-only command is deprecated; "
        "use 'annostat analyze' instead",
        file=sys.stderr,
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
    counts = result["severity_counts"]
    finding_count = sum(counts.values())
    files = write_validation(args.output, result) if finding_count else []
    if not args.quiet:
        status = "PASS" if result["valid"] else "FAIL"
        print(f"Annostat {__version__} | annotation validation {status}")
        if finding_count:
            print(f"  Errors                   {counts['error']:,}")
            print(f"  Warnings                 {counts['warning']:,}")
            print(f"  Findings                 {finding_count:,}")
            print(f"  Output                   {args.output.resolve()}")
            print(f"  Files written            {len(files)}")
        else:
            print("  No validation issues found; no output files were written.")
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
