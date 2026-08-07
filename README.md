# annostat

`annostat` analyzes bacterial genome annotations from matching GFF3 and FASTA files.
It uses only the Python standard library and supports multi-record FASTA files.

## Features

- Counts CDS, RNA types, all feature types, hypothetical proteins, gene annotations,
  and CDS with COG categories
- Exports a compact feature overview as CSV or TSV
- Extracts every CDS on the correct strand as nucleotide and amino-acid multi-FASTA
- Calculates codon usage percentages, start-codon counts, and individual COG category counts
- Saves COG-category and CDS-length plots as dependency-free SVG images

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
- `cog_categories.svg` and `cds_lengths.svg`: visualizations

## Tests

```powershell
python3 -m unittest discover -v
```
