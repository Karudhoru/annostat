"""Comparative annotation profiling for two or more bacterial genomes."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from dataclasses import dataclass
from html import escape
from pathlib import Path
from statistics import median
from typing import Callable, Iterable
from urllib.parse import parse_qs, urlparse

from annostat import __version__
from annostat.analysis import COG_CATEGORY_NAMES, analyze_features
from annostat.models import Feature
from annostat.parsers import parse_fasta, parse_gff
from annostat.plots import write_cog_comparison, write_comparison_overview
from annostat.qc import feature_quality_findings, quality_summary, sequence_quality_findings
from annostat.sequences import (
    declared_genetic_codes,
    iter_cds_sequences,
    recognized_start_codons,
    resolve_genetic_code,
)
from annostat.validation import validate_annotation


@dataclass(frozen=True, slots=True)
class GenomeInput:
    """Identify one labelled FASTA/GFF3 pair in a comparison."""

    label: str
    fasta: Path
    gff: Path
    metadata: dict[str, object] | None = None
    genetic_code: int | None = None


def _sha256(path: Path) -> str:
    """Return the hexadecimal SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _cog_ids(feature: Feature) -> set[str]:
    """Return explicit identifiers such as ``COG0123`` from ``Dbxref``."""

    identifiers = set()
    for cross_reference in feature.attributes.get("Dbxref", "").split(","):
        namespace, separator, value = cross_reference.partition(":")
        if separator and namespace == "COG" and value.startswith("COG") and value[3:].isdigit():
            identifiers.add(value)
    return identifiers


def _species_name(organism_name: str | None) -> str | None:
    """Return a conservative binomial species label from an organism name."""

    if not organism_name:
        return None
    words = organism_name.split()
    if len(words) < 2:
        return None
    lowered = {word.lower().strip(".,;()[]") for word in words}
    if lowered & {"uncultured", "unclassified", "unidentified", "environmental", "metagenome"}:
        return None
    genus_index = 1 if words[0].lower() == "candidatus" else 0
    if len(words) <= genus_index + 1:
        return None
    species_word = words[genus_index + 1].lower().strip(".,;()[]")
    if species_word in {"sp", "bacterium", "archaeon", "cf", "aff"}:
        return None
    word_count = 3 if genus_index else 2
    return " ".join(words[:word_count])


def _assembly_display(accession: object) -> str:
    """Render a safe NCBI assembly link when a versioned accession is known."""

    value = str(accession) if accession else ""
    if value.startswith(("GCF_", "GCA_")):
        url = f"https://www.ncbi.nlm.nih.gov/datasets/genome/{value}/"
        return f'<a href="{escape(url)}">{escape(value)}</a>'
    return escape(value or "-")


def _gff_metadata(path: Path) -> dict[str, object]:
    """Read organism and annotation provenance from leading GFF3 comments."""

    metadata: dict[str, object] = {"input_source": "Local files"}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            text = line.strip()
            if text.startswith("# organism "):
                metadata["organism_name"] = text.removeprefix("# organism ").strip()
            elif text.startswith("# annotated with "):
                metadata["annotation_pipeline"] = text.removeprefix("# annotated with ").strip()
            elif text.startswith("#!genome-build-accession "):
                accession = text.removeprefix("#!genome-build-accession ").strip()
                metadata["assembly_accession"] = accession.removeprefix("NCBI_Assembly:")
            elif text.startswith("#!genome-build "):
                metadata["assembly_name"] = text.removeprefix("#!genome-build ").strip()
            elif text.startswith("#!annotation-source "):
                metadata["annotation_provider"] = text.removeprefix("#!annotation-source ").strip()
            elif text.startswith("##species "):
                taxonomy_reference = text.removeprefix("##species ").strip()
                query_id = parse_qs(urlparse(taxonomy_reference).query).get("id", [])
                if query_id and query_id[0].isdigit():
                    metadata["taxonomy_id"] = int(query_id[0])
                else:
                    trailing_id = taxonomy_reference.rstrip("/").rsplit("/", 1)[-1]
                    if trailing_id.isdigit():
                        metadata["taxonomy_id"] = int(trailing_id)
    return metadata


def _taxonomic_relationship(left: dict[str, object], right: dict[str, object]) -> str:
    """Classify a pair using known species labels without guessing missing data."""

    left_species, right_species = left.get("species_name"), right.get("species_name")
    if not left_species or not right_species:
        return "not determined"
    return "same species" if left_species == right_species else "different species"


def _annotation_profile(dataset: GenomeInput) -> dict[str, object]:
    """Calculate the normalized comparison profile for one input dataset."""

    validation = validate_annotation(dataset.fasta, dataset.gff)
    if not validation["valid"]:
        errors = [
            finding for finding in validation["findings"]
            if finding["severity"] == "error"
        ]
        first = errors[0]
        raise ValueError(
            f"dataset {dataset.label!r} failed validation with {len(errors)} error(s): "
            f"{first['rule_id']}: {first['message']}"
        )
    features = list(parse_gff(dataset.gff))
    genome = parse_fasta(dataset.fasta)
    circular_seqids = frozenset(
        feature.seqid
        for feature in features
        if feature.type == "region"
        and feature.attributes.get("Is_circular", "").lower() == "true"
    )
    declared_codes = declared_genetic_codes(features)
    genetic_code = resolve_genetic_code(features, dataset.genetic_code)
    genetic_code_source = (
        "command_line" if dataset.genetic_code is not None
        else "gff3" if declared_codes
        else "default"
    )
    findings = feature_quality_findings(features, genome, circular_seqids)
    start_counts: Counter[str] = Counter()
    cds_lengths: list[int] = []
    for record in iter_cds_sequences(
        features,
        genome,
        circular_seqids,
        start_counts=start_counts,
        genetic_code=genetic_code,
    ):
        cds_lengths.append(record.length)
        findings.extend(
            sequence_quality_findings(
                record,
                genetic_code,
                sequence_length=len(genome[record.feature.seqid]),
                circular=record.feature.seqid in circular_seqids,
            )
        )

    feature_summary = analyze_features(features)
    quality = quality_summary(features, genome, circular_seqids, findings)
    cds_count = int(feature_summary["cds_count"])
    genome_length = sum(map(len, genome.values()))
    gene_symbols = {
        feature.attributes["gene"].strip()
        for feature in features
        if feature.type == "CDS" and feature.attributes.get("gene", "").strip()
    }
    cog_ids = set().union(
        *(_cog_ids(feature) for feature in features if feature.type == "CDS")
    ) if features else set()
    source_counts = Counter(
        feature.source for feature in features if feature.source and feature.source != "."
    )
    cog_counts = feature_summary["cog_category_counts"]
    cog_data_available = bool(cog_counts)
    cog_total = sum(cog_counts.values())
    start_total = sum(start_counts.values())
    recognized_starts = sum(
        start_counts.get(codon, 0) for codon in recognized_start_codons(genetic_code)
    )
    metadata = {**_gff_metadata(dataset.gff), **(dataset.metadata or {})}
    organism_name = metadata.get("organism_name")
    return {
        "label": dataset.label,
        "input_source": metadata.get("input_source", "Local files"),
        "organism_name": organism_name,
        "species_name": _species_name(str(organism_name)) if organism_name else None,
        "taxonomy_id": metadata.get("taxonomy_id"),
        "assembly_accession": metadata.get("assembly_accession"),
        "assembly_name": metadata.get("assembly_name"),
        "assembly_level": metadata.get("assembly_level"),
        "annotation_provider": metadata.get("annotation_provider"),
        "annotation_release": metadata.get("annotation_release"),
        "annotation_pipeline": metadata.get("annotation_pipeline"),
        "genetic_code": genetic_code,
        "genetic_code_source": genetic_code_source,
        "checkm_completeness": metadata.get("checkm_completeness"),
        "checkm_contamination": metadata.get("checkm_contamination"),
        "inputs": {"fasta": str(dataset.fasta), "gff3": str(dataset.gff)},
        "sha256": {"fasta": _sha256(dataset.fasta), "gff3": _sha256(dataset.gff)},
        "validation": {
            "valid": validation["valid"],
            "ruleset_version": validation["ruleset_version"],
            "scientific_fingerprint": validation["scientific_fingerprint"],
            "severity_counts": validation["severity_counts"],
        },
        "annotation_sources": dict(sorted(source_counts.items())),
        "genome_length": genome_length,
        "sequence_count": len(genome),
        "total_features": int(feature_summary["total_features"]),
        "cds_count": cds_count,
        "rna_count": sum(feature_summary["rna_counts"].values()),
        "gc_percent": quality["genome_gc_percent"],
        "coding_density_percent": quality["coding_density_percent"],
        "cds_per_mb": 1_000_000 * cds_count / genome_length if genome_length else 0,
        "median_cds_length": median(cds_lengths) if cds_lengths else 0,
        "hypothetical_percent": 100 * int(feature_summary["hypothetical_cds_count"]) / cds_count if cds_count else 0,
        "gene_name_percent": 100 * int(feature_summary["cds_with_gene_count"]) / cds_count if cds_count else 0,
        "cog_coverage_percent": (
            100 * int(feature_summary["cds_with_cog_count"]) / cds_count
            if cds_count and cog_data_available else None
        ),
        "cog_data_available": cog_data_available,
        "recognized_start_percent": 100 * recognized_starts / cds_count if cds_count else 0,
        "qc_warning_count": quality["severity_counts"].get("warning", 0),
        "qc_issue_counts": quality["issue_counts"],
        "cog_category_percentages": {
            category: 100 * count / cog_total if cog_total else 0
            for category, count in sorted(cog_counts.items())
        },
        "start_codon_percentages": {
            codon: 100 * count / start_total if start_total else 0
            for codon, count in sorted(start_counts.items())
        },
        "gene_symbols": sorted(gene_symbols),
        "cog_ids": sorted(cog_ids),
    }


def _jaccard(left: Iterable[str], right: Iterable[str]) -> float | None:
    """Return Jaccard similarity, or ``None`` when both sets are empty."""

    left_set, right_set = set(left), set(right)
    union = left_set | right_set
    return len(left_set & right_set) / len(union) if union else None


def _jensen_shannon(left: dict[str, float], right: dict[str, float]) -> float | None:
    """Return base-2 Jensen-Shannon distance for two percentage profiles."""

    labels = set(left) | set(right)
    if not labels or not sum(left.values()) or not sum(right.values()):
        return None
    left_total, right_total = sum(left.values()), sum(right.values())
    p = {label: left.get(label, 0) / left_total for label in labels}
    q = {label: right.get(label, 0) / right_total for label in labels}
    midpoint = {label: (p[label] + q[label]) / 2 for label in labels}

    def divergence(values: dict[str, float]) -> float:
        """Measure one profile's Kullback-Leibler divergence from the midpoint."""

        return sum(
            value * math.log2(value / midpoint[label])
            for label, value in values.items()
            if value
        )

    return math.sqrt((divergence(p) + divergence(q)) / 2)


_NUMERIC_METRICS = (
    "genome_length", "sequence_count", "total_features", "cds_count", "rna_count",
    "gc_percent", "coding_density_percent", "cds_per_mb", "median_cds_length",
    "hypothetical_percent", "gene_name_percent", "cog_coverage_percent",
    "recognized_start_percent", "qc_warning_count",
)


def _pairwise_rows(profiles: list[dict[str, object]]) -> list[dict[str, object]]:
    """Build deterministic pairwise metric and annotation-overlap rows."""

    rows = []
    for left_index, left in enumerate(profiles):
        for right in profiles[left_index + 1:]:
            row: dict[str, object] = {
                "dataset_a": left["label"],
                "dataset_b": right["label"],
                "taxonomic_relationship": _taxonomic_relationship(left, right),
                "gene_symbol_jaccard": _jaccard(left["gene_symbols"], right["gene_symbols"]),
                "cog_id_jaccard": (
                    _jaccard(left["cog_ids"], right["cog_ids"])
                    if left["cog_data_available"] and right["cog_data_available"]
                    else None
                ),
                "cog_profile_distance": _jensen_shannon(
                    left["cog_category_percentages"], right["cog_category_percentages"]
                ),
                "start_profile_distance": _jensen_shannon(
                    left["start_codon_percentages"], right["start_codon_percentages"]
                ),
            }
            for metric in _NUMERIC_METRICS:
                if left[metric] is None or right[metric] is None:
                    row[f"{metric}_delta"] = None
                    row[f"{metric}_percent_delta"] = None
                    continue
                left_value, right_value = float(left[metric]), float(right[metric])
                row[f"{metric}_delta"] = right_value - left_value
                row[f"{metric}_percent_delta"] = (
                    100 * (right_value - left_value) / left_value if left_value else None
                )
            rows.append(row)
    return rows


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    """Write dictionaries as a TSV with stable columns and blank null values."""

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if value is None else value for key, value in row.items()})


def _profile_table_rows(profiles: list[dict[str, object]]) -> list[dict[str, object]]:
    """Select scalar profile values for the public dataset metrics table."""

    return [
        {
            "dataset": profile["label"],
            "organism": profile["organism_name"],
            "species": profile["species_name"],
            "input_source": profile["input_source"],
            "assembly_accession": profile["assembly_accession"],
            "annotation": profile["annotation_provider"] or profile["annotation_pipeline"],
            "genetic_code": profile["genetic_code"],
            "genetic_code_source": profile["genetic_code_source"],
            **{metric: profile[metric] for metric in _NUMERIC_METRICS},
            "cog_data_available": profile["cog_data_available"],
            "cog_coverage_percent": (
                profile["cog_coverage_percent"] if profile["cog_data_available"] else None
            ),
        }
        for profile in profiles
    ]


def _render_report(summary: dict[str, object], plot_paths: list[Path]) -> str:
    """Render a compact, self-contained comparative HTML report."""

    profiles = summary["datasets"]
    species_count = len(summary["species"])
    scope_label = {
        "mixed_species": f"{species_count} identified species",
        "single_species": "one identified species",
        "not_determined": "taxonomic scope not fully determined",
    }[summary["taxonomic_scope"]]
    overview_rows = "".join(
        "<tr>"
        f'<td><strong>{escape(str(profile["label"]))}</strong></td>'
        f'<td>{escape(str(profile["organism_name"] or "Not provided"))}</td>'
        f'<td>{escape(str(profile["input_source"]))}</td>'
        f'<td>{profile["genetic_code"]}</td>'
        f'<td>{profile["genome_length"]:,} bp</td><td>{profile["cds_count"]:,}</td>'
        f'<td>{profile["gc_percent"]:.2f}%</td><td>{profile["coding_density_percent"]:.2f}%</td>'
        f'<td>{_display_cog_coverage(profile)}</td></tr>'
        for profile in profiles
    )
    pair_rows = "".join(
        "<tr>"
        f'<td>{escape(str(row["dataset_a"]))} / {escape(str(row["dataset_b"]))}</td>'
        f'<td><span class="scope">{escape(str(row["taxonomic_relationship"]))}</span></td>'
        f'<td>{_display_similarity(row["gene_symbol_jaccard"])}</td>'
        f'<td>{_display_similarity(row["cog_id_jaccard"])}</td>'
        f'<td>{_display_distance(row["cog_profile_distance"])}</td>'
        f'<td>{_display_distance(row["start_profile_distance"])}</td></tr>'
        for row in summary["pairwise_comparisons"]
    )
    warning_panel = ""
    if summary["warnings"]:
        warning_items = "".join(
            f"<li>{escape(warning)}</li>" for warning in summary["warnings"]
        )
        warning_panel = f'<aside class="notice"><strong>Interpret carefully</strong><ul>{warning_items}</ul></aside>'
    figure_sections = "".join(
        f'<section class="figure">{path.read_text(encoding="utf-8")}</section>'
        for path in plot_paths
    )
    provenance_rows = "".join(
        "<tr>"
        f'<td>{escape(str(profile["label"]))}</td>'
        f'<td>{_assembly_display(profile["assembly_accession"])}</td>'
        f'<td>{escape(str(profile["annotation_provider"] or profile["annotation_pipeline"] or "Not provided"))}</td>'
        f'<td>{escape(Path(profile["inputs"]["fasta"]).name)}</td>'
        f'<td>{escape(Path(profile["inputs"]["gff3"]).name)}</td></tr>'
        for profile in profiles
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Annostat comparison report</title><style>
:root{{--ink:#17231e;--muted:#64746c;--line:#dbe5df;--paper:#fff;--bg:#f4f7f5;--green:#096b4f;--soft:#e7f4ee;--warn:#fff8e6;--warn-line:#e7c45c}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 Inter,Segoe UI,Arial,sans-serif}}
header{{background:#0b533e;color:#fff;padding:20px max(24px,calc((100% - 1180px)/2))}}header span{{float:right;color:#d9eee6}}main{{max-width:1180px;margin:auto;padding:34px 24px 64px}}
h1{{font-size:34px;line-height:1.15;margin:0 0 8px}}h2{{font-size:20px;margin:0 0 14px}}.subtitle{{color:var(--muted);margin:0 0 20px}}
.tags{{display:flex;gap:8px;flex-wrap:wrap;margin-bottom:22px}}.tag,.scope{{background:var(--soft);border-radius:999px;color:var(--green);display:inline-block;font-size:12px;font-weight:700;padding:4px 9px}}.unavailable{{color:var(--muted);font-style:italic}}
.panel,.figure,details,.notice{{background:var(--paper);border:1px solid var(--line);border-radius:12px;margin:16px 0;padding:18px}}.table-wrap{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;white-space:nowrap}}th{{color:var(--muted);font-size:11px;letter-spacing:.04em;text-transform:uppercase}}th,td{{border-bottom:1px solid var(--line);padding:10px;text-align:right}}th:nth-child(-n+3),td:nth-child(-n+3){{text-align:left}}tbody tr:last-child td{{border-bottom:0}}
.figure svg{{display:block;height:auto;max-width:100%;width:100%}}.notice{{background:var(--warn);border-color:var(--warn-line)}}.notice ul{{margin:8px 0 0;padding-left:20px}}.note{{color:var(--muted);font-size:13px}}
details summary{{cursor:pointer;font-weight:700}}details .table-wrap{{margin-top:14px}}.downloads{{display:flex;gap:18px;flex-wrap:wrap}}a{{color:var(--green);font-weight:600}}footer{{border-top:1px solid var(--line);color:var(--muted);font-size:12px;margin-top:28px;padding-top:16px}}
@media(max-width:700px){{header span{{display:block;float:none;margin-top:3px}}main{{padding:26px 14px 48px}}h1{{font-size:28px}}.panel,.figure,details,.notice{{padding:14px}}}}@media print{{body{{background:#fff}}main{{max-width:none;padding:12px}}.panel,.figure,details,.notice{{break-inside:avoid}}}}
</style></head><body><header><strong>Annostat comparison</strong><span>version {escape(str(summary["annostat_version"]))}</span></header><main>
<h1>Genome annotation comparison</h1><p class="subtitle">A normalized comparison of annotated bacterial assemblies. Dataset labels do not imply strain-level relatedness.</p>
<div class="tags"><span class="tag">{len(profiles)} assemblies</span><span class="tag">{escape(scope_label)}</span></div>
{warning_panel}
<section class="panel"><h2>Assemblies</h2><div class="table-wrap"><table><thead><tr><th>Dataset</th><th>Organism</th><th>Input source</th><th>Genetic code</th><th>Genome</th><th>CDS</th><th>GC</th><th>Coding</th><th>COG coverage</th></tr></thead><tbody>{overview_rows}</tbody></table></div></section>
{figure_sections}
<section class="panel"><h2>Pairwise comparison</h2><div class="table-wrap"><table><thead><tr><th>Pair</th><th>Taxonomic scope</th><th>Gene labels</th><th>COG IDs</th><th>COG profile distance</th><th>Start profile distance</th></tr></thead><tbody>{pair_rows}</tbody></table></div>
<p class="note">Jaccard values range from 0 to 1 and describe exact annotation-label overlap. Profile distances range from 0 (same distribution) to 1 (maximally different). Neither measure establishes orthology, ANI, or phylogeny.</p></section>
<details><summary>Input sources and reproducibility</summary><div class="table-wrap"><table><thead><tr><th>Dataset</th><th>Assembly</th><th>Annotation</th><th>FASTA</th><th>GFF3</th></tr></thead><tbody>{provenance_rows}</tbody></table></div><p class="note">Complete paths, file checksums, NCBI metadata, and calculation details are retained in comparison.json.</p></details>
<section class="panel"><h2>Download results</h2><div class="downloads"><a href="tables/dataset_metrics.tsv">Dataset metrics</a><a href="tables/pairwise_comparisons.tsv">Pairwise comparisons</a><a href="comparison.json">Complete JSON</a>{''.join(f'<a href="plots/{escape(path.name)}">{escape(path.stem.replace("_", " ").title())} (SVG)</a>' for path in plot_paths)}</div></section>
<footer>Generated locally by Annostat. The report contains no external scripts or network resources.</footer>
</main></body></html>"""


def _display_similarity(value: float | None) -> str:
    """Format an optional similarity value for HTML."""

    return f"{value:.3f}" if value is not None else "n/a"


def _display_cog_coverage(profile: dict[str, object]) -> str:
    """Format COG coverage without presenting missing annotation as zero."""

    if not profile["cog_data_available"]:
        return '<span class="unavailable">Not available</span>'
    return f'{float(profile["cog_coverage_percent"]):.2f}%'


def _display_distance(value: float | None) -> str:
    """Format an optional profile distance for HTML."""

    return f"{value:.3f}" if value is not None else "insufficient data"


def run_comparison(
    datasets: list[GenomeInput],
    output_dir: Path,
    progress: Callable[[str], None] | None = None,
) -> dict[str, object]:
    """Compare labelled genome annotations and write the report package."""

    if len(datasets) < 2:
        raise ValueError("comparison requires at least two genomes")
    labels = [dataset.label for dataset in datasets]
    if len(labels) != len(set(labels)):
        raise ValueError("comparison genome labels must be unique")
    notify = progress or (lambda message: None)
    notify("Analyzing and quality-checking input genomes")
    profiles = [_annotation_profile(dataset) for dataset in datasets]
    notify("Calculating normalized and pairwise comparisons")
    pairwise = _pairwise_rows(profiles)
    source_signatures = {tuple(profile["annotation_sources"]) for profile in profiles}
    species_values = [profile["species_name"] for profile in profiles]
    known_species = sorted({str(species) for species in species_values if species})
    taxonomic_scope = (
        "mixed_species" if len(known_species) > 1
        else "single_species" if len(known_species) == 1 and all(species_values)
        else "not_determined"
    )
    warnings = []
    if taxonomic_scope == "mixed_species":
        warnings.append(
            "The inputs represent different species. Functional profiles remain comparable, but gene-label overlap must not be interpreted as strain-level similarity."
        )
    if len(source_signatures) > 1:
        warnings.append(
            "Annotation source fields differ between datasets; pipeline differences may explain part of the observed variation."
        )
    if len({profile["genetic_code"] for profile in profiles}) > 1:
        warnings.append(
            "Translation tables differ between datasets. Each genome was translated and quality-checked with its own declared or selected genetic code."
        )
    if any(not profile["cog_data_available"] for profile in profiles):
        warnings.append(
            "At least one annotation source does not provide COG categories. COG coverage, overlap, distances, and comparative plots involving that dataset are not available."
        )
    summary = {
        "annostat_version": __version__,
        "analysis_type": "comparative_annotation_profile",
        "taxonomic_scope": taxonomic_scope,
        "species": known_species,
        "datasets": profiles,
        "pairwise_comparisons": pairwise,
        "warnings": warnings,
    }

    tables_dir, plots_dir = output_dir / "tables", output_dir / "plots"
    tables_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    notify("Writing comparison tables and scientific figures")
    _write_rows(tables_dir / "dataset_metrics.tsv", _profile_table_rows(profiles))
    _write_rows(tables_dir / "pairwise_comparisons.tsv", pairwise)
    overview_plot = plots_dir / "dataset_overview.svg"
    cog_plot = plots_dir / "cog_comparison.svg"
    write_comparison_overview(overview_plot, profiles)
    cog_plot.unlink(missing_ok=True)
    plot_paths = [overview_plot]
    if all(profile["cog_data_available"] for profile in profiles):
        write_cog_comparison(cog_plot, profiles, COG_CATEGORY_NAMES)
        plot_paths.append(cog_plot)
    notify("Building the offline comparative report")
    (output_dir / "comparison.html").write_text(
        _render_report(summary, plot_paths), encoding="utf-8"
    )
    (output_dir / "comparison.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary
