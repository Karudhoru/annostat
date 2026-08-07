# annostat

`annostat` analyzes bacterial genome annotations from matching GFF3 and FASTA files.
It uses only the Python standard library and supports multi-record FASTA files.

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
  and start codons; the COG chart is omitted when the source has no COG data
- Streams CDS processing and FASTA generation to reduce runtime and peak memory
- Compares two or more local or NCBI-hosted annotations using normalized genome,
  annotation-completeness, COG, start-codon, and QC profiles
- Reports pairwise Jensen-Shannon profile distances and exact gene-symbol/COG-ID
  overlap without presenting annotation labels as inferred orthologs
- Downloads versioned RefSeq or GenBank assemblies through the official NCBI
  Datasets command-line tool when network connectivity is requested

## Requirements

- Python 3.10 or newer
- Matching sequence identifiers in the GFF3 and FASTA inputs

## Input assumptions and scope

- GFF3 records must contain nine tab-separated fields with 1-based inclusive
  coordinates. The FASTA identifier is the first whitespace-delimited header token.
- CDS strand and phase are honored. Circular-origin CDS coordinates extending past
  the sequence length are supported when a `region` feature declares
  `Is_circular=true`.
- A CDS split across multiple GFF3 rows is currently analyzed as multiple rows; this
  is appropriate for the supplied bacterial annotations but should be considered
  before using eukaryotic or heavily fragmented gene models.
- COG coverage describes categories actually present in the annotation source. An
  annotation with no COG category fields is reported as unavailable, not as 0%.

## Installation

For development, install annostat in editable mode from the repository root:

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

## Testing

Run the complete standard-library test suite from the repository root with the
same interpreter used to install AnnStat:

```bash
python -m compileall -q annostat tests
python -m unittest discover -v
```

Then run one supplied dataset end to end and open `results/report.html` locally:

```bash
annostat -f data/GCF_000007145.1.fna -g data/GCF_000007145.1.gff3 -o results --table-format tsv
```

## Usage

After installation, run annostat from any directory:

```bash
annostat -f data/GCF_000007145.1.fna -g data/GCF_000007145.1.gff3 -o results
```

Without installing, module execution remains available from the repository root:
`python3 -m annostat --help`.

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
assemblies or any mixture of local and NCBI inputs. AnnStat downloads genome
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
not supply COG fields, AnnStat reports COG coverage and pairwise COG statistics
as `Not available` and omits the COG comparison plot; missing data is never
presented as biological zero coverage.

For users who only need the source files, the standalone downloader remains
available:

```bash
annostat fetch GCF_000007145.1 GCF_001050915.2 --output ncbi_data
```

The `NCBI_API_KEY` environment variable is honored by the official client. Local
analysis never requires a network connection, and AnnStat does not transmit
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

## Example run

```bash
annostat \
  -f data/GCF_000007145.1.fna \
  -g data/GCF_000007145.1.gff3 \
  -o results \
  --table-format tsv
```

Example output:

```text
annostat 0.5.0 | bacterial genome annotation analysis
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
  ATG/GTG/TTG starts                     4,294 (99.98%)
  QC review                1 warning, 872 informational findings
--------------------------------------------------------
  Output                   /path/to/annostat/results
  Files written            12
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
