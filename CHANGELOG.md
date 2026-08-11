# Changelog

All notable changes to Annostat are documented in this file. The Python package
and command-line executable retain the lowercase technical name `annostat`.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0]

### Added

- Separate `annostat validate` and `annostat inspect` workflows so deterministic
  structural failures are not conflated with context-dependent biological review
  findings.
- Versioned validation rules with stable identifiers, scientific rationale,
  source URLs, canonical JSON/TSV output, input SHA-256 hashes, and reproducible
  scientific fingerprints.
- `annostat summarize` for MultiQC-style aggregation of completed inspections
  into deterministic cohort HTML, TSV, and JSON reports.
- Correct biological joining and phase validation for multipart CDS features on
  both strands.
- Translation tables 4, 11, and 25 with annotation-aware selection from GFF3
  `transl_table` attributes and safe table-11 fallback.
- Position-specific `transl_except` support for selenocysteine, pyrrolysine,
  termination, and conservative `OTHER`-to-`X` translation.
- A top-level command overview for `annostat --help` while retaining the legacy
  `annostat -f ... -g ...` inspection form.
- Per-genome `--genetic-code LABEL TABLE` overrides for comparative analysis.

### Changed

- Completed the command-specific help for validation, cohort summarization, and
  NCBI fetching with examples, option descriptions, and exit-status guidance.
- Replaced assignment-specific wording in active user and API documentation with
  tool-focused language.
- Declared the MIT license, supported Python versions, scientific audience, and
  repository links in the installable package metadata.
- Inspection and comparison now reject explicit genetic-code selections that
  conflict with a translation table declared by the annotation.
- Start-codon metrics use every initiator recognized by the selected NCBI table,
  rather than counting only ATG, GTG, and TTG.
- Comparison reports retain the selected genetic code and its provenance for
  every dataset.

### Fixed

- Mistyped subcommands now produce a direct unknown-command error with the valid
  choices instead of being interpreted as an incomplete legacy inspection.
- Alternative initiator codons are converted to methionine only for CDS records
  with a complete 5-prime boundary; partial N-termini preserve the ordinary codon
  translation.
- Terminal `transl_except` residues are retained after normal stop-codon trimming.
- Translation-table selection now reads NCBI `region` declarations as well as CDS
  declarations.
- Anonymous CDS rows remain distinct in feature and CDS-density counts, while
  repeated explicit IDs are still joined as multipart features.
- Partial-feature QC suppresses only the affected biological boundary check.
- Comparative analysis now performs structural validation before calculating any
  scientific profiles.
- Reusing an output directory removes obsolete Annostat table/filter variants
  without deleting unrelated user files.
- Circular features may extend their end into virtual origin-spanning coordinates,
  but their start must remain within the real sequence coordinate range.
- Multipart CDS median lengths in comparative reports now use the complete
  joined biological feature rather than the first GFF3 segment.
- NCBI assembly accessions prefixed with `NCBI_Assembly:` and taxonomy URLs using
  `wwwtax.cgi?id=` are normalized for report metadata and links.
- Table-25 annotations no longer generate false internal-stop warnings or
  incorrect proteins when no manual CLI option is supplied.

## [0.5.0]

### Added

- `annostat compare` for normalized, QC-aware comparison of two or more labelled
  FASTA/GFF3 annotation datasets.
- Dataset metrics and pairwise TSV tables containing scalar deltas,
  Jensen-Shannon COG/start-profile distances, and exact gene-symbol/COG-ID
  Jaccard overlap.
- A self-contained comparative HTML report with SHA-256 input provenance,
  interpretation warnings, and explicit scientific limitations.
- Two focused comparison figures: a normalized annotation overview and an
  adaptive two-genome COG difference chart or multi-genome heatmap.
- Optional NCBI connectivity through `annostat fetch` and repeatable
  `annostat compare --reference`, using the official NCBI Datasets CLI.
- Versioned GCF/GCA accession validation and path-safe extraction of downloaded
  NCBI data packages.

### Changed

- Updated the package version to 0.5.0.
- Expanded CLI help and README documentation for local and accession-based
  comparative workflows.
- Dataset labels are taxonomically neutral; local GFF3 and NCBI assembly
  metadata now identify organisms and same-species or different-species pairs.
- Simplified the comparative HTML report, moved detailed reproducibility data
  into a collapsed section and JSON, and corrected wide-table/figure layout.
- External NCBI assemblies are now presented as first-class comparison
  references rather than only as standalone downloads.

### Fixed

- Missing COG annotations from sources such as current NCBI RefSeq GFF3
  packages are reported as unavailable rather than as zero coverage.
- COG overlap, distance, and comparison plots are omitted whenever an involved
  dataset lacks COG data, preventing misleading comparisons against zero.

### Scientific scope

- Gene-symbol overlap is explicitly reported as annotation-label similarity,
  not orthology or core/accessory genome inference.
- Annostat does not implement approximate ANI, phylogeny, or homolog clustering;
  those analyses remain delegated to specialized tools.

## [0.4.0]

### Added

- A consolidated `tables/annotation_issues.csv` containing conservative
  annotation-quality findings with severity, coordinates, related features, and
  human-readable details.
- Detection of overlapping and fully contained CDS pairs, explicitly marked
  pseudogenes, partial features, adjacent duplicate annotations, non-triplet
  CDS lengths, internal stop codons, and ambiguous coding bases.
- Structural-RNA completeness checks for 5S, 16S, and 23S rRNAs and tRNAs for
  the 20 standard amino acids.
- Genome GC percentage, non-duplicated coding density, CDS density, and
  severity-grouped QC counts in `summary.json`, the CLI, and the HTML report.
- Optional filtered CDS exports through `--min-cds-length`,
  `--max-cds-length`, `--require-cog`, and `--exclude-hypothetical`.

### Changed

- Filtered records are written separately under `filtered/`; the complete
  assignment-required analysis and outputs always remain unchanged.
- The HTML report now includes annotation-quality findings and documents active
  filter criteria without adding another default plot.
- Updated the package version to 0.4.0.

## [0.3.0]

### Added

- A self-contained offline `report.html` containing the main annotation metrics,
  provenance, performance information, codon summary, and embedded SVG figures.
- `--profile` for peak Python-memory measurement and per-stage terminal timings.
- Structured `tables/`, `sequences/`, and `plots/` output directories.
- Mean and median markers in the CDS-length distribution.
- Counts and percentages in the COG and start-codon charts.
- The five most-used codons in `summary.json` and the HTML report.

### Changed

- CDS extraction, translation, codon accumulation, and FASTA export now operate
  as a streaming pipeline.
- Translation and codon counting share one codon traversal.
- FASTA records are written in buffered blocks rather than one line per write.
- Feature and COG statistics are calculated in a single pass.
- Start codons are grouped as ATG, GTG, TTG, and Other in the visualization while
  the complete raw counts remain available in `start_codons.csv`.
- Updated the package version to 0.3.0.

### Removed

- The default 64-codon heatmap. Codon usage remains fully available in the
  assignment-required `codon_usage.csv` percentage table.

### Performance

- The supplied `GCF_000007145.1` dataset averaged 0.480 seconds versus 0.562
  seconds for version 0.2.0 in the same local benchmark, a 14.6% improvement
  while also generating the new HTML report.

### Compatibility

- All six table and FASTA outputs shared with version 0.2.0 remain byte-for-byte
  identical for both supplied datasets.

## [0.2.0]

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

## [0.1.0]

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
