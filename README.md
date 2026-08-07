# annostat

`annostat` analyzes bacterial genome annotations from matching GFF3 and FASTA files.
It uses only the Python standard library and supports multi-record FASTA files.

## Features

- Counts CDS, RNA types, all feature types, hypothetical proteins, gene annotations,
  and CDS with COG categories
- Exports a compact feature overview as CSV or TSV
- Extracts every CDS on the correct strand as nucleotide and amino-acid multi-FASTA
- Calculates codon usage percentages, start-codon counts, and individual COG category counts
- Saves four dependency-free SVG charts for COG categories, CDS lengths,
  start codons, and codon usage

## Requirements

- Python 3.10 or newer
- Matching sequence identifiers in the GFF3 and FASTA inputs

## Usage

Run directly from a repository checkout:

```powershell
python3 -m annostat -f data/GCF_000007145.1.fna -g data/GCF_000007145.1.gff3 -o results
```

Direct execution is also supported: `python3 annostat/cli.py --help`.

Use `--table-format tsv` for a tab-separated feature overview. Install the package
with `python -m pip install .` to make the equivalent `annostat` command available.

The output directory contains:

- `summary.json`: required feature and annotation metrics
- `features.csv` or `features.tsv`: compact feature table
- `cds_nucleotide.fasta` and `cds_protein.fasta`: extracted CDS sequences
- `codon_usage.csv`, `start_codons.csv`, and `cog_categories.csv`: analysis tables
- `cog_categories.svg` and `cds_lengths.svg`: annotation visualizations
- `start_codons.svg` and `codon_usage.svg`: codon visualizations

## Example run

```bash
python3 annostat/cli.py \
  -f data/GCF_000007145.1.fna \
  -g data/GCF_000007145.1.gff3 \
  -o results \
  --table-format tsv
```

Example output:

```text
annostat 0.2.0 | bacterial genome annotation analysis
FASTA: /path/to/annostat/data/GCF_000007145.1.fna
GFF3:  /path/to/annostat/data/GCF_000007145.1.gff3

[1/5] Reading GFF3 annotations and FASTA sequences
[2/5] Extracting and translating CDS sequences
[3/5] Calculating feature, COG, and codon statistics
[4/5] Writing tables, summaries, and FASTA files
[5/5] Rendering scientific visualizations

Analysis summary
--------------------------------------------------------
  Genome size                              5,076,188 bp
  Sequences                              1 (1 circular)
  Features                                        4,498
  CDS                                             4,295
  RNA features                                      181
  Hypothetical CDS                         619 (14.41%)
  COG-annotated CDS                      2,881 (67.08%)
  ATG/GTG/TTG starts                     4,294 (99.98%)
--------------------------------------------------------
  Output                   /path/to/annostat/results
  Files written            11
  Completed in             1.16 seconds
```

Runtime depends on the input size and system. The supplied bacterial genome
typically completes in a few seconds.
