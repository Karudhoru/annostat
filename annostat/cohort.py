"""Cohort aggregation for completed Annostat inspection results."""

from __future__ import annotations

import csv
import json
from collections import Counter
from html import escape
from pathlib import Path
from annostat import __version__


COHORT_SCHEMA_VERSION = "1.0"


def discover_summaries(inputs: list[Path]) -> list[Path]:
    """Return unique ``summary.json`` files found beneath explicit inputs."""

    discovered: set[Path] = set()
    for raw_path in inputs:
        path = Path(raw_path)
        if path.is_file():
            if path.name != "summary.json":
                raise ValueError(f"expected an Annostat summary.json file, got {path}")
            discovered.add(path.resolve())
        elif path.is_dir():
            discovered.update(candidate.resolve() for candidate in path.rglob("summary.json"))
        else:
            raise ValueError(f"input path does not exist: {path}")
    if not discovered:
        raise ValueError("no Annostat summary.json files found")
    return sorted(discovered, key=lambda item: item.as_posix().casefold())


def _percent(part: int, whole: int) -> float | None:
    """Return a percentage or ``None`` when the denominator is unavailable."""

    return 100 * part / whole if whole else None


def _sample_row(path: Path) -> dict[str, object]:
    """Read and normalize one single-genome summary."""

    try:
        summary = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read Annostat summary {path}: {error}") from error
    required = ("genome_length", "cds_count", "rna_counts", "quality_control")
    missing = [key for key in required if key not in summary]
    if missing:
        raise ValueError(f"{path} is missing required keys: {', '.join(missing)}")

    cds_count = int(summary["cds_count"])
    hypothetical = int(summary.get("hypothetical_cds_count", 0))
    cog_count = int(summary.get("cds_with_cog_count", 0))
    quality = summary["quality_control"]
    validation = summary.get("validation", {})
    cog_available = bool(summary.get("cog_data_available", False))
    return {
        "sample": path.parent.name,
        "summary_path": str(path),
        "annostat_version": str(summary.get("annostat_version", "unknown")),
        "scientific_fingerprint": str(summary.get("scientific_fingerprint", "")),
        "genome_length_bp": int(summary["genome_length"]),
        "sequence_count": len(summary.get("sequence_ids", [])),
        "cds_count": cds_count,
        "rna_count": sum(int(value) for value in summary["rna_counts"].values()),
        "genome_gc_percent": float(quality["genome_gc_percent"]),
        "coding_density_percent": float(quality["coding_density_percent"]),
        "hypothetical_cds_percent": _percent(hypothetical, cds_count),
        "cog_coverage_percent": _percent(cog_count, cds_count) if cog_available else None,
        "inspection_warning_count": int(quality.get("severity_counts", {}).get("warning", 0)),
        "validation_error_count": int(validation.get("severity_counts", {}).get("error", 0)),
        "validation_warning_count": int(validation.get("severity_counts", {}).get("warning", 0)),
        "validation_status": (
            "pass" if validation.get("valid") is True
            else "fail" if validation.get("valid") is False
            else "not_available"
        ),
        "input_sha256": validation.get("input_sha256", {}),
    }


def build_cohort(inputs: list[Path]) -> dict[str, object]:
    """Aggregate normalized metrics from one or more Annostat outputs."""

    rows = [_sample_row(path) for path in discover_summaries(inputs)]
    labels = [str(row["sample"]) for row in rows]
    label_counts = Counter(labels)
    duplicate_labels = sorted(label for label, count in label_counts.items() if count > 1)
    if duplicate_labels:
        raise ValueError(
            "sample directory names must be unique; duplicates: " + ", ".join(duplicate_labels)
        )
    rows.sort(key=lambda row: str(row["sample"]).casefold())
    versions = sorted({str(row["annostat_version"]) for row in rows})
    return {
        "schema_version": COHORT_SCHEMA_VERSION,
        "annostat_version": __version__,
        "sample_count": len(rows),
        "source_annostat_versions": versions,
        "mixed_source_versions": len(versions) > 1,
        "samples": rows,
    }


def _display(value: object, digits: int = 2) -> str:
    """Format optional numeric values for TSV and HTML output."""

    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def write_cohort_table(path: Path, cohort: dict[str, object]) -> None:
    """Write the normalized sample matrix as deterministic TSV."""

    fields = [
        "sample", "validation_status", "genome_length_bp", "sequence_count",
        "cds_count", "rna_count", "genome_gc_percent", "coding_density_percent",
        "hypothetical_cds_percent", "cog_coverage_percent",
        "validation_error_count", "validation_warning_count",
        "inspection_warning_count", "annostat_version", "summary_path",
        "scientific_fingerprint",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in cohort["samples"]:
            writer.writerow({key: _display(row.get(key)) for key in fields})


def _status_class(status: str) -> str:
    """Return a conservative visual class for a validation status."""

    return {"pass": "pass", "fail": "fail"}.get(status, "unknown")


def render_cohort_html(cohort: dict[str, object]) -> str:
    """Render a self-contained, script-free cohort report."""

    rows = []
    for sample in cohort["samples"]:
        status = str(sample["validation_status"])
        cells = (
            escape(str(sample["sample"])),
            f'<span class="status {_status_class(status)}">{escape(status.replace("_", " "))}</span>',
            f'{int(sample["genome_length_bp"]):,}',
            f'{int(sample["cds_count"]):,}',
            _display(sample["genome_gc_percent"]),
            _display(sample["coding_density_percent"]),
            _display(sample["hypothetical_cds_percent"]),
            _display(sample["cog_coverage_percent"]),
            str(sample["validation_error_count"]),
            str(sample["inspection_warning_count"]),
        )
        rows.append("<tr>" + "".join(f"<td>{cell}</td>" for cell in cells) + "</tr>")
    version_note = (
        '<aside class="notice">Inputs were produced by multiple Annostat versions; '
        "compare metrics with the recorded version metadata.</aside>"
        if cohort["mixed_source_versions"] else ""
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Annostat cohort report</title><style>
:root{{--ink:#17231e;--muted:#66756d;--line:#dce6e0;--paper:#fff;--bg:#f3f7f5;--green:#087f5b;--red:#b42318}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font-family:Segoe UI,Arial,sans-serif;line-height:1.45}}
header{{background:#0c563f;color:#fff;padding:22px max(28px,calc((100% - 1320px)/2));font-size:22px;font-weight:700}}
main{{max-width:1320px;margin:auto;padding:38px 28px 70px}}h1{{font-size:36px;margin:0 0 4px}}.subtitle{{color:var(--muted);margin:0 0 24px}}
.panel,.notice{{background:var(--paper);border:1px solid var(--line);border-radius:14px;padding:20px;overflow-x:auto}}.notice{{margin-bottom:18px;border-left:5px solid #d99a00}}
table{{width:100%;border-collapse:collapse;font-size:14px}}th,td{{padding:10px 9px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}}th{{color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.04em}}th:first-child,td:first-child{{text-align:left;font-weight:600}}
.status{{font-weight:700;text-transform:uppercase;font-size:11px}}.pass{{color:var(--green)}}.fail{{color:var(--red)}}.unknown{{color:var(--muted)}}
footer{{color:var(--muted);font-size:12px;margin-top:24px}}@media print{{body{{background:#fff}}main{{max-width:none;padding:15px}}}}
</style></head><body><header>Annostat</header><main><h1>Cohort annotation QC</h1>
<p class="subtitle">{int(cohort["sample_count"]):,} inspection result(s), ordered by sample identifier.</p>{version_note}
<section class="panel"><table><thead><tr><th>Sample</th><th>Validation</th><th>Genome bp</th><th>CDS</th><th>GC %</th><th>Coding %</th><th>Hypothetical %</th><th>COG %</th><th>Errors</th><th>Inspection warnings</th></tr></thead><tbody>{''.join(rows)}</tbody></table></section>
<footer>NA means the source annotation did not provide that metric. Biological inspection warnings are review candidates, not automatic errors. Generated locally by Annostat {escape(__version__)}.</footer>
</main></body></html>"""


def write_cohort(output_dir: Path, cohort: dict[str, object]) -> list[Path]:
    """Write the complete deterministic cohort result package."""

    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "cohort.json"
    table_path = output_dir / "cohort.tsv"
    html_path = output_dir / "cohort.html"
    json_path.write_text(
        json.dumps(cohort, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_cohort_table(table_path, cohort)
    html_path.write_text(render_cohort_html(cohort), encoding="utf-8")
    return [json_path, table_path, html_path]
