"""Offline HTML reporting for annostat analyses."""

from __future__ import annotations

from html import escape
from pathlib import Path


def _percentage(part: int, whole: int) -> str:
    """Format a part-to-whole ratio as a percentage, including zero totals."""

    return f"{100 * part / whole:.2f}%" if whole else "0.00%"


def _metric(label: str, value: str, note: str) -> str:
    """Render one escaped summary metric card."""

    return (
        '<article class="metric">'
        f'<div class="metric-label">{escape(label)}</div>'
        f'<div class="metric-value">{escape(value)}</div>'
        f'<div class="metric-note">{escape(note)}</div>'
        "</article>"
    )


def _embedded_svg(path: Path, heading: str) -> str:
    """Read an SVG plot and wrap it in a report figure section."""

    svg = path.read_text(encoding="utf-8")
    return f'<section class="figure"><h2>{escape(heading)}</h2>{svg}</section>'


def render_html_report(
    summary: dict[str, object],
    plot_paths: list[tuple[Path, str]],
) -> str:
    """Render a self-contained HTML analysis report with inline SVG figures.

    User-controlled paths and labels are HTML-escaped, and the generated report
    has no external scripts, fonts, stylesheets, or network dependencies.
    """

    cds_count = int(summary["cds_count"])
    rna_count = sum(summary["rna_counts"].values())
    hypothetical = int(summary["hypothetical_cds_count"])
    cog_annotated = int(summary["cds_with_cog_count"])
    start_counts = summary["start_codon_counts"]
    standard_starts = sum(start_counts.get(codon, 0) for codon in ("ATG", "GTG", "TTG"))
    input_files = summary["input_files"]
    performance = summary["performance"]
    stage_rows = "".join(
        f"<tr><td>{escape(label.replace('_', ' ').title())}</td><td>{seconds:.4f} s</td></tr>"
        for label, seconds in performance["stage_seconds"].items()
    )
    peak_memory = performance.get("peak_memory_bytes")
    peak_memory_text = f"{peak_memory / 1024 / 1024:.2f} MiB" if peak_memory else "Run with --profile to measure"
    top_codons = "".join(
        f'<li><strong>{escape(item["codon"])}</strong><span>{item["count"]:,} ({item["percentage"]:.2f}%)</span></li>'
        for item in summary["top_codons"]
    )
    metrics = "".join(
        (
            _metric("Genome size", f'{int(summary["genome_length"]):,} bp', f'{len(summary["sequence_ids"]):,} sequence records'),
            _metric("Features", f'{int(summary["total_features"]):,}', f'{cds_count:,} coding sequences'),
            _metric("COG coverage", _percentage(cog_annotated, cds_count), f'{cog_annotated:,} annotated CDS'),
            _metric("Recognized starts", _percentage(standard_starts, cds_count), f'{standard_starts:,} ATG/GTG/TTG'),
        )
    )
    figures = "".join(_embedded_svg(path, heading) for path, heading in plot_paths)
    generated_files = [
        path for path in summary["output_files"] if path != "report.html"
    ]
    file_links = "".join(
        f'<li><a href="{escape(path)}">{escape(path)}</a></li>' for path in generated_files
    )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>annostat bacterial annotation report</title>
<style>
:root{{--ink:#17231e;--muted:#66756d;--line:#dce6e0;--paper:#fff;--background:#f3f7f5;--green:#087f5b;--green-soft:#dff5ec}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--background);color:var(--ink);font-family:Inter,Segoe UI,Arial,sans-serif;line-height:1.5}}
header{{background:#0c563f;color:#fff;padding:22px max(28px,calc((100% - 1240px)/2))}} header strong{{font-size:22px}} header span{{float:right;color:#d8f4e8}}
main{{max-width:1240px;margin:auto;padding:42px 28px 70px}} h1{{font-size:38px;line-height:1.15;margin:0 0 8px}} h2{{font-size:21px;margin:0 0 16px}} .subtitle{{color:var(--muted);margin-bottom:28px}}
.metrics{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:22px}} .metric,.panel,.figure{{background:var(--paper);border:1px solid var(--line);border-radius:14px;padding:20px}}
.metric-label{{color:var(--muted);font-size:12px;font-weight:700;letter-spacing:.05em;text-transform:uppercase}} .metric-value{{font-size:27px;font-weight:700;margin:7px 0 2px}} .metric-note{{color:var(--muted);font-size:13px}}
.summary-grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin:22px 0}} dl{{display:grid;grid-template-columns:190px 1fr;gap:9px 16px;margin:0}} dt{{color:var(--muted)}} dd{{margin:0;font-weight:600;overflow-wrap:anywhere}}
.codons{{list-style:none;padding:0;margin:0;display:grid;gap:9px}} .codons li{{display:flex;justify-content:space-between;border-bottom:1px solid var(--line);padding-bottom:7px}}
.figures{{display:grid;gap:20px}} .figure svg{{width:100%;height:auto}} .figure h2{{margin-bottom:8px}}
table{{width:100%;border-collapse:collapse}} td{{border-bottom:1px solid var(--line);padding:8px 4px}} td:last-child{{text-align:right;font-variant-numeric:tabular-nums}}
.files{{columns:2;margin:0;padding-left:20px}} a{{color:var(--green)}} footer{{color:var(--muted);font-size:12px;margin-top:32px;border-top:1px solid var(--line);padding-top:18px}}
@media(max-width:800px){{.metrics{{grid-template-columns:1fr 1fr}}.summary-grid{{grid-template-columns:1fr}}}} @media(max-width:520px){{.metrics{{grid-template-columns:1fr}}dl{{grid-template-columns:1fr}}header span{{float:none;display:block}}.files{{columns:1}}}}
@media print{{body{{background:#fff}}main{{max-width:none;padding:20px}}.metric,.panel,.figure{{break-inside:avoid}}}}
</style>
</head>
<body>
<header><strong>annostat</strong><span>version {escape(str(summary["annostat_version"]))}</span></header>
<main>
<h1>Bacterial annotation report</h1>
<p class="subtitle">{escape(Path(input_files["gff3"]).name)} analyzed against {escape(Path(input_files["fasta"]).name)}</p>
<section class="metrics">{metrics}</section>
<section class="summary-grid">
  <article class="panel"><h2>Annotation summary</h2><dl>
    <dt>RNA features</dt><dd>{rna_count:,}</dd>
    <dt>Hypothetical CDS</dt><dd>{hypothetical:,} ({_percentage(hypothetical, cds_count)})</dd>
    <dt>Gene names</dt><dd>{int(summary["cds_with_gene_count"]):,} ({_percentage(int(summary["cds_with_gene_count"]), cds_count)})</dd>
    <dt>Complete codons</dt><dd>{int(summary["complete_codon_count"]):,}</dd>
    <dt>Circular sequences</dt><dd>{len(summary["circular_sequence_ids"]):,}</dd>
  </dl></article>
  <article class="panel"><h2>Most-used codons</h2><ol class="codons">{top_codons}</ol></article>
</section>
<section class="figures">{figures}</section>
<section class="summary-grid">
  <article class="panel"><h2>Performance</h2><table>{stage_rows}<tr><td>Total measured stages</td><td>{performance["total_seconds"]:.4f} s</td></tr><tr><td>Peak Python memory</td><td>{escape(peak_memory_text)}</td></tr></table></article>
  <article class="panel"><h2>Methods and provenance</h2><dl>
    <dt>Translation table</dt><dd>NCBI bacterial genetic code 11</dd>
    <dt>FASTA</dt><dd>{escape(str(input_files["fasta"]))}</dd>
    <dt>GFF3</dt><dd>{escape(str(input_files["gff3"]))}</dd>
    <dt>Codon calculation</dt><dd>Pooled across complete CDS codons</dd>
  </dl></article>
</section>
<section class="panel"><h2>Generated data files</h2><ul class="files">{file_links}</ul></section>
<footer>Generated locally by annostat {escape(str(summary["annostat_version"]))}. The report contains no external scripts, fonts, or network resources.</footer>
</main>
</body>
</html>
"""
