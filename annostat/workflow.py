"""Reusable single-annotation analysis workflow."""

from __future__ import annotations

import hashlib
import json
import tracemalloc
from collections import Counter
from pathlib import Path
from time import perf_counter
from typing import Callable

from annostat import __version__
from annostat.analysis import COG_CATEGORY_NAMES, analyze_features
from annostat.filtering import CdsFilter
from annostat.output import (
    write_annotation_findings,
    write_cds_fastas,
    write_codon_usage,
    write_count_table,
    write_overview,
    write_summary,
)
from annostat.plots import write_bar_chart, write_histogram
from annostat.qc import (
    feature_quality_findings,
    quality_summary,
    sequence_quality_findings,
)
from annostat.report import render_html_report
from annostat.sequences import (
    declared_genetic_codes,
    iter_cds_sequences,
    recognized_start_codons,
    resolve_genetic_code,
)
from annostat.validation import load_and_validate_annotation, write_validation


ANALYSIS_SCHEMA_VERSION = "1.0"


def _scientific_fingerprint(summary: dict[str, object]) -> str:
    """Hash scientific outputs while excluding paths and performance timings."""

    payload = {
        key: value
        for key, value in summary.items()
        if key not in {
            "annostat_version",
            "input_files",
            "output_files",
            "performance",
            "schema_version",
            "scientific_fingerprint",
        }
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


def _run_analysis(
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

    notify("Reading GFF3 annotations and FASTA sequences")
    stage_started = perf_counter()
    validated = load_and_validate_annotation(fasta_path, gff_path)
    validation = validated.validation
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
    features = validated.features
    genome = validated.genome
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
    summary["schema_version"] = ANALYSIS_SCHEMA_VERSION
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
    validation_files = write_validation(output_dir / "validation", validation)
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
        "No complete first codon": len(cds_lengths) - sum(start_counts.values()),
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
        *(path.relative_to(output_dir).as_posix() for path in validation_files),
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
    """Run the complete annotation-analysis workflow and return its summary."""

    if profile:
        tracemalloc.start()
    try:
        return _run_analysis(
            fasta_path,
            gff_path,
            output_dir,
            table_format,
            progress,
            profile,
            cds_filter,
            genetic_code,
        )
    finally:
        if profile:
            tracemalloc.stop()
