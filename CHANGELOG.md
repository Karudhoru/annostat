# Changelog

All notable changes to `annostat` are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-08-07

### Added

- `--version` for reporting the installed annostat version.
- `--quiet` for scripts and automated workflows that do not need console output.
- Five-stage progress reporting for parsing, extraction, analysis, export, and plotting.
- A structured terminal summary showing genome size, feature counts, annotation
  percentages, standard start-codon usage, output location, and runtime.
- Version and input-file provenance in `summary.json`.
- A start-codon usage chart.
- A 64-codon usage heatmap grouped by codon position.
- Full descriptions for COG functional-category labels.
- Accessible titles, descriptions, roles, and data tooltips in generated SVG files.
- Regression tests for version output, quiet mode, progress reporting, summaries,
  and SVG validity.

### Changed

- Redesigned the COG-category chart with descriptive labels, gridlines, and a
  consistent scientific-report theme.
- Replaced the CDS-length bar representation with a true vertical histogram.
- Expanded CLI help with stable command naming, clearer metavariables, defaults,
  examples, and output descriptions.
- Expanded the README with a complete example command and representative output.
- Increased the generated visualization set from two charts to four.

### Compatibility

- Core summaries and all six scientific table/FASTA outputs remain identical to
  version 0.1.0 for both supplied reference datasets.
- The additional reporting and visualizations introduce no measurable runtime
  regression in the project benchmark.

## [0.1.0] - 2026-08-07

### Added

- Modular GFF3 and multi-record FASTA parsers.
- Structured feature and CDS sequence models.
- Counts for CDS, RNA types, feature types, hypothetical proteins, gene
  annotations, and COG-annotated CDS.
- Individual handling of multi-letter COG categories such as `COG:OE`.
- Strand-aware CDS extraction using 1-based inclusive GFF3 coordinates.
- Circular chromosome and plasmid origin-crossing sequence extraction.
- Bacterial genetic-code translation with alternative start-codon handling.
- Nucleotide and amino-acid multi-FASTA exports.
- CSV or TSV feature overviews with empty values preserved.
- Codon-usage percentages, start-codon counts, and COG-category tables.
- COG-category and CDS-length SVG visualizations.
- An `argparse` command-line interface and installable `annostat` entry point.
- Unit and end-to-end tests covering parsing, coordinates, translation, analysis,
  circular sequences, direct CLI execution, and generated outputs.

### Fixed

- Invalid GFF3 fields now produce contextual validation errors.
- Empty FASTA headers, empty records, duplicate identifiers, and sequence data
  before a header are rejected with clear errors.
- Direct execution through `python3 annostat/cli.py` resolves package imports
  correctly.
