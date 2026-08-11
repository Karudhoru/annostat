# Annostat

**Annostat** analyzes bacterial genome annotations from matching GFF3 and FASTA
files. Local analysis uses only the Python standard library and supports multi-record FASTA files.

## Features

- Counts CDS, RNA types, all feature types, hypothetical proteins, gene annotations,
  and CDS with COG categories
- Exports a compact feature overview as CSV or TSV
- Extracts every CDS on the correct strand as nucleotide and amino-acid multi-FASTA
- Calculates codon usage percentages, start-codon counts, and individual COG category counts
- Reports GC content, coding density, CDS overlaps and containment, pseudogenes,
  partial features, structural-RNA completeness, frame issues, and ambiguous bases
- Supports optional CDS exports filtered by length, COG annotation, and
  hypothetical-protein status without changing the complete analysis
- Builds a self-contained offline HTML report with provenance and performance timings
- Saves three focused, dependency-free SVG charts for COG categories, CDS lengths,
  and table-aware start codons; the COG chart is omitted when the source has no
  COG data
- Streams CDS processing and FASTA generation to reduce runtime and peak memory
- Compares two or more local or NCBI-hosted annotations using normalized genome,
  annotation-completeness, COG, start-codon, and QC profiles
- Reports pairwise Jensen-Shannon profile distances and exact gene-symbol/COG-ID
  overlap without presenting annotation labels as inferred orthologs
- Downloads versioned RefSeq or GenBank assemblies through the official NCBI
  Datasets command-line tool when network connectivity is requested
- Separates deterministic FASTA/GFF3 validation from biological inspection, with
  stable rule IDs, severity levels, schema/ruleset versions, and source citations
- Records SHA-256 hashes for reproducibility and writes canonical JSON plus TSV
  validation results suitable for CI pipelines
- Correctly joins repeated-ID multipart CDS features in biological order on both
  strands and validates phase continuity before translation
- Supports NCBI prokaryotic translation tables 4, 11, and 25, automatically
  selecting a consistent CDS or `region` `transl_table` declaration when present
- Applies position-specific NCBI `transl_except` annotations for selenocysteine,
  pyrrolysine, explicitly annotated termination codons, and `OTHER` residues
  represented conservatively as `X`
- Aggregates completed inspections into a script-free cohort HTML report,
  normalized TSV matrix, and machine-readable JSON

## Required and optional components

| Component | Status | Purpose |
|---|---|---|
| Python 3.10 or newer | Required | Runs Annostat and its test suite |
| Genome FASTA (`.fna`, `.fa`, or `.fasta`) | Required | Supplies the genome sequences |
| Matching GFF3 annotation (`.gff` or `.gff3`) | Required | Supplies feature coordinates and annotations |
| Matching FASTA/GFF3 sequence identifiers | Required | Connects every annotated feature to its sequence |
| NCBI Datasets CLI | Optional | Downloads GCF/GCA assemblies for `fetch` or `compare --reference` |
| Internet connection | Optional | Needed only when downloading NCBI data |

No third-party Python package is required for local analysis, report generation,
tables, FASTA exports, or SVG plots.

## Assignment-required functionality

Annostat implements the required deliverables from the exercise specification:

- Counts CDS, RNA features by type, hypothetical CDS, CDS with gene names, CDS
  with COG categories, and other feature types.
- Writes a CSV or TSV feature overview containing ID, type, start, stop, gene,
  product, and strand, with missing values left empty.
- Exports all CDS as nucleotide and translated amino-acid multi-FASTA files.
- Calculates pooled codon usage percentages, start-codon counts, and detailed
  COG-category statistics, including multi-category assignments.
- Generates at least two meaningful scientific plots; a normal run produces CDS
  length, start-codon, and—when available—COG-category SVG figures.
- Provides an `argparse` command-line interface and modular Python structure with
  documented functions, classes, and methods.

Filtering, quality control, HTML reports, performance profiling, comparative
analysis, and NCBI connectivity are additional features beyond the core requirements.

## Input assumptions and scope

- GFF3 records must contain nine tab-separated fields with 1-based inclusive
  coordinates. The FASTA identifier is the first whitespace-delimited header token.
- CDS strand and phase are honored. Circular-origin CDS coordinates extending past
  the sequence length are supported when a `region` feature declares
  `Is_circular=true`.
- Repeated CDS IDs are interpreted as one multipart biological feature. Segments
  must use the same sequence ID, strand, and Parent; phase continuity is checked,
  and segments are joined in biological order before translation.
- COG coverage describes categories actually present in the annotation source. An
  annotation with no COG category fields is reported as unavailable, not as 0%.
- Annostat uses the consistent CDS or `region` `transl_table` declared in GFF3.
  When the annotation declares no table, it falls back to bacterial table 11. `--genetic-code` is an
  explicit override for annotations without that metadata; an override that
  conflicts with a declared table is rejected instead of silently mistranslating
  CDS features. Annostat never guesses a genetic code from sequence composition.

## Installation

For development, install Annostat in editable mode from the repository root:

```bash
python -m pip install -e .
annostat --version
```

Editable mode keeps the `annostat` command connected to this checkout, so source
changes take effect without reinstalling. If Bash has cached command locations,
run `hash -r` once after installation.

For a regular installation that copies the current version into the active Python
environment, use:

```bash
python -m pip install .
```

## Step-by-step guide

### 1. Open the repository

Run all setup commands from the directory containing `pyproject.toml`:

```bash
cd /path/to/exam_3
```

### 2. Select one Python installation

Use the same interpreter for installation, testing, and execution. On Linux or
WSL this is commonly `python3`; on Windows it may be `python` or `py -3`:

```bash
python3 --version
```

The reported version must be Python 3.10 or newer.

### 3. Create and activate a virtual environment

Linux or WSL:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Windows PowerShell:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
```

### 4. Install Annostat

Editable installation is recommended while developing:

```bash
python -m pip install -e .
annostat --version
annostat --help
```

The expected version banner begins with `Annostat`. If the `annostat` command is
not found, use `python -m annostat` with the same arguments.

### 5. Check the inputs

Confirm that you have one genome FASTA and its matching GFF3 annotation. Sequence
IDs must agree—for example, a FASTA record named `chromosome` must be referenced
as `chromosome` in the first GFF3 column.

### 6. Validate the paired files

Run deterministic integrity checks before biological interpretation:

```bash
annostat validate \
  -f data/GCF_000007145.1.fna \
  -g data/GCF_000007145.1.gff3 \
  -o validation
```

The command writes `validation.json` and `validation.tsv`. Exit status `0` means
no structural errors were found; status `1` means the configured threshold was
reached. Use `--fail-on warning` for strict CI, or `--fail-on never` when the
report should never stop a pipeline. A validation pass means the paired files are
internally interpretable; it does not prove that every biological annotation is
correct.

### 7. Inspect one annotation

```bash
annostat inspect \
  -f data/GCF_000007145.1.fna \
  -g data/GCF_000007145.1.gff3 \
  -o results
```

Here `-f` selects the required FASTA, `-g` selects the required GFF3, and `-o`
chooses the output directory. Add `--table-format tsv` when a TSV overview is
preferred over CSV. The original command without the `inspect` word remains a
backward-compatible alias. Translation-table selection is automatic when the
GFF3 CDS or `region` rows contain `transl_table`; otherwise table 11 is used. For a local
annotation without that attribute, select a supported table explicitly with
`--genetic-code 4`, `--genetic-code 11`, or `--genetic-code 25`.

### 8. Review the result package

Open `results/report.html` in a browser first. Then use:

- `summary.json` for the complete machine-readable summary;
- `tables/` for feature, codon, COG, start-codon, and QC data;
- `sequences/` for nucleotide and protein CDS FASTA files; and
- `plots/` for editable SVG figures.

### 9. Summarize a cohort

After inspecting several genomes into separate result directories, aggregate the
entire parent directory:

```bash
annostat summarize results/batch -o cohort
```

The output contains `cohort.html`, `cohort.tsv`, and `cohort.json`. Samples are
ordered deterministically. Missing source-dependent values such as COG coverage
remain `NA`/`null`, never biological zero. Multiple source Annostat versions are
retained and flagged in the report.

### 10. Use optional features

The following sections show filtered exports, comparative analysis, NCBI
references, performance profiling, and complete output descriptions.

Without installing, module execution remains available from the repository root:
`python -m annostat --help`.

Use `--table-format tsv` for a tab-separated feature overview. Add `--profile`
to measure peak Python memory and print per-stage timings.

### Filtered CDS exports

Filters create an additional `filtered/` directory while retaining every normal
analysis file. Criteria are combined, so a CDS must satisfy all enabled options:

```bash
annostat \
  -f genome.fna \
  -g annotations.gff3 \
  -o results \
  --min-cds-length 300 \
  --require-cog \
  --exclude-hypothetical
```

Available criteria are `--min-cds-length`, `--max-cds-length`, `--require-cog`,
and `--exclude-hypothetical`.

### Comparative analysis

Compare two or more labelled local genomes with repeatable `--genome` options:

```bash
annostat compare \
  --genome genome-a data/GCF_000007145.1.fna data/GCF_000007145.1.gff3 \
  --genome genome-b data/GCF_001050915.2.fna data/GCF_001050915.2.gff3 \
  --output comparison
```

Each genome is structurally validated and independently uses its declared
`transl_table` or the table-11 fallback. For an annotation without
`transl_table`, a per-genome override uses
the same label supplied to `--genome` (or the reference accession):

```bash
annostat compare \
  --genome sr1 sr1.fna sr1.gff3 \
  --genome ecoli ecoli.fna ecoli.gff3 \
  --genetic-code sr1 25 \
  --output comparison
```

The comparison package is intentionally compact:

```text
comparison/
|-- comparison.html
|-- comparison.json
|-- tables/
|   |-- dataset_metrics.tsv
|   `-- pairwise_comparisons.tsv
`-- plots/
    |-- dataset_overview.svg
    `-- cog_comparison.svg (only when every input supplies COG categories)
```

`dataset_metrics.tsv` provides raw genome/feature counts together with GC,
coding density, CDS per Mb, annotation coverage, start-codon usage, and QC
warnings. `pairwise_comparisons.tsv` contains metric deltas, percentage deltas,
Jensen-Shannon distances for COG and start-codon profiles, and Jaccard overlap
for exact gene symbols and explicit COG identifiers.

For two datasets, the COG figure shows percentage-point differences. For three
or more datasets, it becomes a dataset-by-category heatmap. Exact gene-symbol
overlap describes annotation labels only; it is not orthology, ANI, phylogeny,
or a core/accessory pangenome analysis.

### Compare against external NCBI data

NCBI access is optional and uses NCBI's maintained `datasets` executable. Follow
the official [installation instructions](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/),
then add an external reference directly to a local comparison:

```bash
annostat compare \
  --genome local-sample sample.fna sample.gff3 \
  --reference GCF_000007145.1 \
  --output comparison
```

`--reference` is repeatable, so a comparison may contain only external
assemblies or any mixture of local and NCBI inputs. Annostat downloads genome
FASTA, GFF3, and assembly metadata under `comparison/external_inputs/`, then
runs every input through the same comparison engine:

```bash
annostat compare \
  --reference GCF_000007145.1 \
  --reference GCF_001050915.2 \
  --output comparison
```

The report identifies organisms and whether each pair belongs to the same or
different species when the input metadata makes that determination possible.
It never assumes that labelled datasets are strains. NCBI assembly names,
taxonomy identifiers, annotation providers, and release metadata are retained
in `comparison.json`. Versioned GCF/GCA accessions in the HTML provenance table
link to their NCBI Datasets assembly pages.

COG results are source-dependent. Bakta annotations commonly include COG
categories, while current NCBI RefSeq GFF3 packages may not. When a source does
not supply COG fields, Annostat reports COG coverage and pairwise COG statistics
as `Not available` and omits the COG comparison plot; missing data is never
presented as biological zero coverage.

For users who only need the source files, the standalone downloader remains
available:

```bash
annostat fetch GCF_000007145.1 GCF_001050915.2 --output ncbi_data
```

The `NCBI_API_KEY` environment variable is honored by the official client. Local
analysis never requires a network connection, and Annostat does not transmit
local FASTA or GFF3 files.

The output directory contains:

```text
results/
├── report.html
├── summary.json
├── tables/
│   ├── features.csv (or features.tsv)
│   ├── codon_usage.csv
│   ├── start_codons.csv
│   ├── cog_categories.csv
│   └── annotation_issues.csv
├── validation/
│   ├── validation.json
│   └── validation.tsv
├── sequences/
│   ├── cds_nucleotide.fasta
│   └── cds_protein.fasta
└── plots/
    ├── cog_categories.svg (when COG categories are available)
    ├── cds_lengths.svg
    └── start_codons.svg
```

When filters are active, `filtered/` additionally contains the selected feature
table and paired nucleotide/protein FASTA files.

`report.html` works without a server or internet connection and embeds the SVG
figures directly. The separate SVG files remain available for editing and use in
documents. Codon usage remains available as the required percentage table; a
dense 64-codon heatmap is intentionally not generated by default.

`annotation_issues.csv` is a review aid rather than an automatic pass/fail
verdict. Short CDS overlaps are common in compact bacterial genomes, so ordinary
overlaps and adjacent duplicate labels are informational. Containment, incomplete
frames, missing structural RNA classes, internal stops, and ambiguous bases are
reported as warnings that merit inspection in their biological context.

## Validation versus inspection

`annostat validate` answers whether the files can be interpreted consistently:
GFF3 syntax, FASTA symbols, matching sequence IDs, coordinate bounds, Parent/ID
relationships, multipart-feature consistency, CDS phase continuity, and valid,
consistent CDS translation-table declarations. These checks are deterministic
and may safely control a CI exit status.

`annostat inspect` first requires structural validation to pass, then reports
biological review candidates such as internal stops, missing terminal codons,
ambiguous CDS bases, overlaps, partial features, pseudogenes, duplicated adjacent
labels, and structural-RNA coverage. These findings are context-dependent and do
not automatically make an annotation invalid. Boundary-specific partial markers
suppress only the affected start- or stop-codon assumption. Explicit `pseudo`,
`transl_except`, and `exception` annotations suppress checks whose assumptions do
not apply.

Every validation result embeds its rule definitions, rationale, source URL,
ruleset version, schema version, Annostat version, and input SHA-256 hashes. JSON
and TSV findings use canonical ordering and contain no timestamp, so repeated runs
on identical paths and bytes produce byte-identical validation artifacts. A
`scientific_fingerprint` in both validation and inspection JSON excludes file
paths and performance timings, allowing scientific results to be compared across
machines and output directories.

## Testing and verification

Run compilation and all automated tests with the interpreter used to install
Annostat:

```bash
python -m compileall -q annostat tests
python -m unittest discover -v
```

For a final end-to-end check, analyze one supplied genome and confirm that
`results/report.html`, `results/summary.json`, the tables, FASTA exports, and SVG
plots are created:

```bash
annostat inspect \
  -f data/GCF_000007145.1.fna \
  -g data/GCF_000007145.1.gff3 \
  -o results \
  --table-format tsv
```

## Example run

```bash
annostat inspect \
  -f data/GCF_000007145.1.fna \
  -g data/GCF_000007145.1.gff3 \
  -o results \
  --table-format tsv
```

Example output:

```text
Annostat 1.0.0 | bacterial genome annotation analysis
FASTA: /path/to/annostat/data/GCF_000007145.1.fna
GFF3:  /path/to/annostat/data/GCF_000007145.1.gff3

[1/5] Reading GFF3 annotations and FASTA sequences
[2/5] Streaming, translating, and writing CDS sequences
[3/5] Calculating statistics and writing analysis tables
[4/5] Rendering scientific visualizations
[5/5] Building the offline HTML report

Analysis summary
--------------------------------------------------------
  Genome size                              5,076,188 bp
  Sequences                              1 (1 circular)
  Features                                        4,498
  CDS                                             4,295
  RNA features                                      181
  Genome GC                                      65.07%
  Coding density                                 85.33%
  Hypothetical CDS                         619 (14.41%)
  COG-annotated CDS                      2,881 (67.08%)
  Recognized starts                      4,294 (99.98%)
  QC review                1 warning, 872 informational findings
--------------------------------------------------------
  Output                   /path/to/annostat/results
  Files written            14
  Completed in             0.56 seconds
```

Runtime depends on the input size and system. The supplied bacterial genome
typically completes in under a second on a modern desktop. Use `--profile` when
you need measured stage timings and peak-memory information; profiling itself
adds overhead and should not be used for speed benchmarks.

## Scientific basis

The QC checks follow the internal-consistency emphasis of the
[NCBI Prokaryotic Genome Annotation Standards](https://www.ncbi.nlm.nih.gov/refseq/annotation_prok/standards/)
and the attribute conventions in the
[NCBI GFF3 documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/).
FASTA validation accepts the NCBI/IUPAC nucleotide alphabet, and translation
tables and initiator codons follow the
[NCBI Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
The report treats overlaps as review candidates because overlapping genes are
common in compact bacterial genomes. Translation uses the bacterial code and
alternative-start behavior described by
[Prodigal](https://pmc.ncbi.nlm.nih.gov/articles/PMC2848648/). The project remains
an annotation-analysis tool rather than a replacement for comprehensive
annotation pipelines such as
[Bakta](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743544/) or
[Prokka](https://pubmed.ncbi.nlm.nih.gov/24642063/).

Comparative profiles are normalized to avoid treating larger genomes as
automatically enriched. Full ortholog and pangenome inference remains the domain
of tools such as [Panaroo](https://pmc.ncbi.nlm.nih.gov/articles/PMC7376924/),
while sequence-level average nucleotide identity should be calculated by a
specialized implementation such as
[FastANI](https://pmc.ncbi.nlm.nih.gov/articles/PMC6269478/).
