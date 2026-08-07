from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from pathlib import Path
from unittest.mock import patch

from annostat.analysis import analyze_codons, analyze_features, cog_categories
from annostat.cli import run_analysis
from annostat.comparison import GenomeInput, _species_name, run_comparison
from annostat.filtering import CdsFilter
from annostat.models import Feature
from annostat.ncbi import _safe_extract, fetch_genomes
from annostat.parsers import parse_attributes, parse_fasta, parse_gff
from annostat.qc import (
    feature_quality_findings,
    quality_summary,
    sequence_quality_findings,
)
from annostat.report import render_html_report
from annostat.sequences import (
    extract_cds_sequences,
    iter_cds_sequences,
    reverse_complement,
    translate_dna,
)


class AnnostatTests(unittest.TestCase):
    def test_attributes_and_multi_letter_cog_categories(self) -> None:
        attributes = parse_attributes(
            "ID=cds1;product=A%20protein;Dbxref=COG:COG1234,COG:OE,GO:1"
        )
        feature = Feature("chr", "test", "CDS", 1, 9, None, "+", 0, attributes)

        self.assertEqual(attributes["product"], "A protein")
        self.assertEqual(cog_categories(feature), ("O", "E"))

    def test_coordinates_strand_phase_and_translation(self) -> None:
        plus = Feature("chr", "test", "CDS", 2, 10, None, "+", 0, {"ID": "plus"})
        minus = Feature("chr", "test", "CDS", 2, 10, None, "-", 0, {"ID": "minus"})
        records = extract_cds_sequences([plus, minus], {"chr": "AATGAAATAGC"})

        self.assertEqual(records[0].nucleotide, "ATGAAATAG")
        self.assertEqual(records[0].protein, "MK")
        self.assertEqual(records[1].nucleotide, reverse_complement("ATGAAATAG"))
        self.assertEqual(translate_dna("GTGAAATAA"), "MK")

    def test_cds_can_wrap_across_a_circular_sequence_origin(self) -> None:
        feature = Feature("plasmid", "test", "CDS", 8, 13, None, "+", 0, {"ID": "wrap"})

        record = extract_cds_sequences(
            [feature], {"plasmid": "AAACCCATG"}, frozenset({"plasmid"})
        )[0]

        self.assertEqual(record.nucleotide, "TGAAAC")

    def test_circular_origin_overlap_and_density_are_normalized(self) -> None:
        features = [
            Feature("plasmid", "test", "CDS", 8, 13, None, "+", 0, {"ID": "wrap"}),
            Feature("plasmid", "test", "CDS", 3, 5, None, "+", 0, {"ID": "origin"}),
        ]
        genome = {"plasmid": "AAACCCATG"}
        circular = frozenset({"plasmid"})

        findings = feature_quality_findings(features, genome, circular)
        summary = quality_summary(features, genome, circular, findings)

        overlaps = [finding for finding in findings if finding.issue_type == "cds_overlap"]
        self.assertEqual(len(overlaps), 1)
        self.assertIn("2 bp overlap", overlaps[0].details)
        self.assertAlmostEqual(summary["coding_density_percent"], 7 / 9 * 100)

    def test_required_counts_and_codon_usage(self) -> None:
        features = [
            Feature(
                "chr", "test", "CDS", 1, 9, None, "+", 0,
                {"ID": "one", "gene": "abc", "product": "hypothetical protein", "Dbxref": "COG:OE"},
            ),
            Feature("chr", "test", "tRNA", 10, 20, None, "+", None, {"ID": "rna"}),
        ]
        summary = analyze_features(features)
        records = extract_cds_sequences(features, {"chr": "ATGAAATAA" + "A" * 11})
        codons, starts = analyze_codons(records)

        self.assertEqual(summary["cds_count"], 1)
        self.assertEqual(summary["rna_counts"], {"tRNA": 1})
        self.assertEqual(summary["hypothetical_cds_count"], 1)
        self.assertEqual(summary["cog_category_counts"], {"E": 1, "O": 1})
        self.assertEqual(codons, {"ATG": 1, "AAA": 1, "TAA": 1})
        self.assertEqual(starts, {"ATG": 1})

    def test_streaming_translation_accumulates_identical_codon_statistics(self) -> None:
        features = [
            Feature("chr", "test", "CDS", 1, 9, None, "+", 0, {"ID": "one"}),
            Feature("chr", "test", "CDS", 10, 18, None, "+", 0, {"ID": "two"}),
        ]
        genome = {"chr": "ATGAAATAAGTGCCCTAG"}
        expected_records = extract_cds_sequences(features, genome)
        expected_codons, expected_starts = analyze_codons(expected_records)
        codons: Counter[str] = Counter()
        starts: Counter[str] = Counter()

        streamed_records = list(
            iter_cds_sequences(
                features, genome, codon_counts=codons, start_counts=starts
            )
        )

        self.assertEqual(streamed_records, expected_records)
        self.assertEqual(codons, expected_codons)
        self.assertEqual(starts, expected_starts)

    def test_feature_analysis_consumes_a_single_pass_iterable(self) -> None:
        feature = Feature("chr", "test", "CDS", 1, 9, None, "+", 0, {"ID": "one"})
        features = iter((feature,))

        summary = analyze_features(features)

        self.assertEqual(summary["total_features"], 1)
        self.assertEqual(summary["cds_count"], 1)

    def test_cds_filter_combines_all_enabled_criteria(self) -> None:
        cds_filter = CdsFilter(
            min_length=9,
            max_length=12,
            require_cog=True,
            exclude_hypothetical=True,
        )
        matching = Feature(
            "chr", "test", "CDS", 1, 9, None, "+", 0,
            {"ID": "match", "product": "enzyme", "Dbxref": "COG:J"},
        )
        hypothetical = Feature(
            "chr", "test", "CDS", 10, 18, None, "+", 0,
            {"ID": "hyp", "product": "hypothetical protein", "Dbxref": "COG:J"},
        )
        missing_cog = Feature(
            "chr", "test", "CDS", 19, 27, None, "+", 0,
            {"ID": "no-cog", "product": "enzyme"},
        )

        self.assertTrue(cds_filter.active)
        self.assertTrue(cds_filter.matches(matching))
        self.assertFalse(cds_filter.matches(hypothetical))
        self.assertFalse(cds_filter.matches(missing_cog))

    def test_annotation_quality_checks_are_conservative_and_structured(self) -> None:
        features = [
            Feature(
                "chr", "test", "CDS", 1, 12, None, "+", 0,
                {"ID": "outer", "gene": "abc", "partial": "true"},
            ),
            Feature(
                "chr", "test", "CDS", 4, 9, None, "-", 0,
                {"ID": "inner", "gene": "def"},
            ),
            Feature(
                "other", "test", "CDS", 1, 9, None, "+", 0,
                {"ID": "duplicate-one", "gene": "dup"},
            ),
            Feature(
                "other", "test", "CDS", 20, 28, None, "+", 0,
                {"ID": "duplicate-two", "gene": "dup"},
            ),
            Feature(
                "chr", "test", "tRNA", 13, 15, None, "+", None,
                {"ID": "pseudo", "pseudo": "True", "product": "tRNA-Ala"},
            ),
        ]
        features.extend(
            Feature(
                "rna", "test", "rRNA", index * 10 + 1, index * 10 + 9,
                None, "+", None, {"ID": kind, "product": f"{kind} ribosomal RNA"},
            )
            for index, kind in enumerate(("5S", "16S", "23S"))
        )
        features.extend(
            Feature(
                "rna", "test", "tRNA", 40 + index * 3, 42 + index * 3,
                None, "+", None, {"ID": amino_acid, "product": f"tRNA-{amino_acid}"},
            )
            for index, amino_acid in enumerate(
                (
                    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
                    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
                )
            )
        )

        findings = feature_quality_findings(features)
        issue_counts = Counter(finding.issue_type for finding in findings)
        self.assertEqual(issue_counts["contained_cds"], 1)
        self.assertEqual(issue_counts["adjacent_duplicate_annotation"], 1)
        self.assertEqual(issue_counts["pseudogene"], 1)
        self.assertEqual(issue_counts["partial_feature"], 1)
        self.assertNotIn("missing_rrna_type", issue_counts)
        self.assertNotIn("missing_trna_amino_acid", issue_counts)

        record = extract_cds_sequences(
            [Feature("problem", "test", "CDS", 1, 10, None, "+", 0, {"ID": "bad"})],
            {"problem": "ATGTAANNNN"},
        )[0]
        sequence_findings = sequence_quality_findings(record)
        sequence_types = {finding.issue_type for finding in sequence_findings}
        self.assertEqual(
            sequence_types,
            {"non_triplet_cds", "internal_stop_codon", "ambiguous_cds_bases"},
        )

        summary = quality_summary(
            features,
            {"chr": "A" * 30, "other": "GC" * 20, "rna": "A" * 120},
            frozenset(),
            findings + sequence_findings,
        )
        self.assertGreater(summary["coding_density_percent"], 0)
        self.assertGreater(summary["genome_gc_percent"], 0)
        self.assertEqual(summary["finding_count"], len(findings) + len(sequence_findings))

    def test_end_to_end_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fasta = root / "genome.fna"
            gff = root / "genes.gff3"
            output = root / "results"
            fasta.write_text(">chr description\nATGAAATAA\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "chr\ttest\tCDS\t1\t9\t.\t+\t0\tID=cds1;product=test;Dbxref=COG:J\n",
                encoding="utf-8",
            )

            summary = run_analysis(
                fasta,
                gff,
                output,
                "tsv",
                profile=True,
                cds_filter=CdsFilter(min_length=9, require_cog=True),
            )

            self.assertEqual(summary["cds_count"], 1)
            expected = {
                "report.html", "summary.json", "tables/features.tsv",
                "tables/codon_usage.csv", "tables/start_codons.csv",
                "tables/cog_categories.csv", "sequences/cds_nucleotide.fasta",
                "sequences/cds_protein.fasta", "plots/cog_categories.svg",
                "plots/cds_lengths.svg", "plots/start_codons.svg",
                "tables/annotation_issues.csv", "filtered/features.tsv",
                "filtered/cds_nucleotide.fasta", "filtered/cds_protein.fasta",
            }
            self.assertEqual(
                {path.relative_to(output).as_posix() for path in output.rglob("*") if path.is_file()},
                expected,
            )
            with (output / "tables" / "features.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["ID"], "cds1")
            with (output / "filtered" / "features.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                filtered_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["ID"] for row in filtered_rows], ["cds1"])
            self.assertIn(
                ">cds1 test",
                (output / "filtered" / "cds_nucleotide.fasta").read_text(),
            )
            written_summary = json.loads((output / "summary.json").read_text())
            self.assertEqual(written_summary["cds_count"], 1)
            self.assertRegex(written_summary["annostat_version"], r"^\d+\.\d+\.\d+$")
            self.assertEqual(written_summary["input_files"]["fasta"], str(fasta))
            self.assertEqual(written_summary["top_codons"][0]["codon"], "ATG")
            self.assertEqual(len(written_summary["output_files"]), 15)
            self.assertEqual(written_summary["filtered_export"]["selected_cds_count"], 1)
            self.assertIn("coding_density_percent", written_summary["quality_control"])
            self.assertIn("stage_seconds", written_summary["performance"])
            self.assertIn("report_generation", written_summary["performance"]["stage_seconds"])
            self.assertGreater(written_summary["performance"]["peak_memory_bytes"], 0)
            report = (output / "report.html").read_text(encoding="utf-8")
            self.assertIn("Bacterial annotation report", report)
            self.assertIn("Most-used codons", report)
            self.assertIn("Annotation quality findings", report)
            self.assertIn("Filtered CDS export", report)
            self.assertIn("Report Generation", report)
            self.assertNotIn("<script", report)
            self.assertNotIn("https://", report)
            self.assertEqual(report.count("<svg "), 3)
            self.assertNotIn("Codon usage heatmap", report)
            for plot_name in ("cog_categories.svg", "cds_lengths.svg", "start_codons.svg"):
                root_element = ET.parse(output / "plots" / plot_name).getroot()
                self.assertTrue(root_element.tag.endswith("svg"))
                self.assertIsNotNone(root_element.find("{http://www.w3.org/2000/svg}title"))

            written_summary["input_files"] = {
                "fasta": "<script>alert(1)</script>.fna",
                "gff3": "unsafe&annotation.gff3",
            }
            escaped_report = render_html_report(
                written_summary,
                [
                    (output / "plots" / "cog_categories.svg", "COG categories"),
                    (output / "plots" / "cds_lengths.svg", "CDS lengths"),
                    (output / "plots" / "start_codons.svg", "Start codons"),
                ],
            )
            self.assertNotIn("<script>alert(1)</script>", escaped_report)
            self.assertIn("unsafe&amp;annotation.gff3", escaped_report)

    def test_single_report_treats_missing_cog_annotations_as_unavailable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fasta = root / "genome.fna"
            gff = root / "genes.gff3"
            output = root / "results"
            fasta.write_text(">chr\nATGAAATAA\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\nchr\ttest\tCDS\t1\t9\t.\t+\t0\tID=cds1;product=test\n",
                encoding="utf-8",
            )

            summary = run_analysis(fasta, gff, output, "csv")

            self.assertFalse(summary["cog_data_available"])
            self.assertNotIn("plots/cog_categories.svg", summary["output_files"])
            self.assertFalse((output / "plots" / "cog_categories.svg").exists())
            report = (output / "report.html").read_text(encoding="utf-8")
            self.assertIn("Not available", report)
            self.assertEqual(report.count("<svg "), 2)

    def test_comparative_analysis_outputs_normalized_profiles(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            datasets = []
            for label, source, sequence, cog in (
                ("alpha", "pipeline-a", "ATGAAATAA", "J"),
                ("beta", "pipeline-b", "GTGCCCTAA", "E"),
            ):
                fasta = root / f"{label}.fna"
                gff = root / f"{label}.gff3"
                fasta.write_text(f">chr\n{sequence}\n", encoding="utf-8")
                gff.write_text(
                    f"##gff-version 3\n# organism {label.title()} species\n"
                    f"chr\t{source}\tCDS\t1\t9\t.\t+\t0\t"
                    f"ID={label}-cds;gene=shared;product=test;"
                    f"Dbxref=COG:{cog},COG:COG0001\n",
                    encoding="utf-8",
                )
                metadata = {"assembly_accession": "GCF_000000001.1"} if label == "alpha" else None
                datasets.append(GenomeInput(label, fasta, gff, metadata))

            output = root / "comparison"
            summary = run_comparison(datasets, output)

            self.assertEqual(len(summary["datasets"]), 2)
            self.assertEqual(summary["pairwise_comparisons"][0]["gene_symbol_jaccard"], 1)
            self.assertEqual(summary["pairwise_comparisons"][0]["cog_id_jaccard"], 1)
            self.assertEqual(summary["pairwise_comparisons"][0]["cog_profile_distance"], 1)
            self.assertEqual(summary["taxonomic_scope"], "mixed_species")
            self.assertEqual(
                summary["pairwise_comparisons"][0]["taxonomic_relationship"],
                "different species",
            )
            self.assertEqual(len(summary["warnings"]), 2)
            expected = {
                "comparison.html", "comparison.json", "tables/dataset_metrics.tsv",
                "tables/pairwise_comparisons.tsv", "plots/dataset_overview.svg",
                "plots/cog_comparison.svg",
            }
            self.assertEqual(
                {path.relative_to(output).as_posix() for path in output.rglob("*") if path.is_file()},
                expected,
            )
            with (output / "tables" / "dataset_metrics.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["dataset"] for row in rows], ["alpha", "beta"])
            report = (output / "comparison.html").read_text(encoding="utf-8")
            self.assertIn("Genome annotation comparison", report)
            self.assertIn("Alpha species", report)
            self.assertIn("different species", report)
            self.assertIn("https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000000001.1/", report)
            self.assertIn("Dataset Overview (SVG)", report)
            self.assertIn("exact annotation-label overlap", report)
            self.assertNotIn("Methods and provenance", report)
            self.assertNotIn("SHA-256 ", report)
            self.assertEqual(report.count("<svg "), 2)
            for plot_name in ("dataset_overview.svg", "cog_comparison.svg"):
                ET.parse(output / "plots" / plot_name)

    def test_species_names_are_not_inferred_from_placeholder_taxa(self) -> None:
        self.assertIsNone(_species_name("uncultured bacterium"))
        self.assertIsNone(_species_name("Escherichia sp. ABC"))
        self.assertIsNone(_species_name("environmental samples"))
        self.assertEqual(_species_name("Candidatus Liberibacter asiaticus"), "Candidatus Liberibacter asiaticus")

    def test_comparison_validates_input_count_and_unique_labels(self) -> None:
        missing = GenomeInput("same", Path("missing.fna"), Path("missing.gff3"))
        with self.assertRaisesRegex(ValueError, "at least two"):
            run_comparison([missing], Path("unused"))
        with self.assertRaisesRegex(ValueError, "unique"):
            run_comparison([missing, missing], Path("unused"))

    def test_comparison_treats_missing_cog_annotations_as_unavailable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            inputs = []
            for label, dbxref in (("annotated", ";Dbxref=COG:J"), ("external", "")):
                fasta = root / f"{label}.fna"
                gff = root / f"{label}.gff3"
                fasta.write_text(">chr\nATGAAATAA\n", encoding="utf-8")
                gff.write_text(
                    "##gff-version 3\n"
                    f"chr\ttest\tCDS\t1\t9\t.\t+\t0\tID={label};product=test{dbxref}\n",
                    encoding="utf-8",
                )
                inputs.append(GenomeInput(label, fasta, gff))

            output = root / "comparison"
            summary = run_comparison(inputs, output)

            self.assertTrue(summary["datasets"][0]["cog_data_available"])
            self.assertFalse(summary["datasets"][1]["cog_data_available"])
            self.assertIsNone(summary["datasets"][1]["cog_coverage_percent"])
            pair = summary["pairwise_comparisons"][0]
            self.assertIsNone(pair["cog_id_jaccard"])
            self.assertIsNone(pair["cog_profile_distance"])
            self.assertFalse((output / "plots" / "cog_comparison.svg").exists())
            report = (output / "comparison.html").read_text(encoding="utf-8")
            self.assertIn("Not available", report)
            self.assertEqual(report.count("<svg "), 1)
            with (output / "tables" / "dataset_metrics.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[1]["cog_data_available"], "False")
            self.assertEqual(rows[1]["cog_coverage_percent"], "")

    def test_ncbi_fetch_rejects_invalid_accessions_before_networking(self) -> None:
        with self.assertRaisesRegex(ValueError, "versioned GCF_ or GCA_"):
            fetch_genomes(["GCF_123"], Path("unused"))

    def test_ncbi_fetch_uses_official_dataset_package_layout(self) -> None:
        accession = "GCF_000000001.1"
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)

            def fake_download(command, **kwargs):
                archive = Path(command[command.index("--filename") + 1])
                with zipfile.ZipFile(archive, "w") as package:
                    base = f"ncbi_dataset/data/{accession}"
                    package.writestr(f"{base}/{accession}_genomic.fna", ">chr\nATG\n")
                    package.writestr(f"{base}/genomic.gff", "##gff-version 3\n")
                    package.writestr("ncbi_dataset/data/dataset_catalog.json", "{}")
                    package.writestr(
                        "ncbi_dataset/data/assembly_data_report.jsonl",
                        json.dumps(
                            {
                                "accession": accession,
                                "organism": {"organismName": "Example bacterium", "taxId": 42},
                                "assemblyInfo": {"assemblyName": "ASM1", "assemblyLevel": "Complete Genome"},
                                "annotationInfo": {"provider": "NCBI RefSeq"},
                            }
                        ) + "\n",
                    )
                return subprocess.CompletedProcess(command, 0, "", "")

            with (
                patch("annostat.ncbi.shutil.which", return_value="datasets"),
                patch("annostat.ncbi.subprocess.run", side_effect=fake_download) as run,
            ):
                fetched = fetch_genomes([accession], root / "download")

            self.assertEqual(fetched[0].accession, accession)
            self.assertTrue(fetched[0].fasta.is_file())
            self.assertTrue(fetched[0].gff.is_file())
            self.assertEqual(fetched[0].metadata["organism_name"], "Example bacterium")
            self.assertEqual(fetched[0].metadata["annotation_provider"], "NCBI RefSeq")
            command = run.call_args.args[0]
            self.assertEqual(command[:5], ["datasets", "download", "genome", "accession", accession])
            self.assertIn("genome,gff3", command)

    def test_ncbi_zip_extraction_rejects_path_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            archive = root / "unsafe.zip"
            with zipfile.ZipFile(archive, "w") as package:
                package.writestr("../escape.txt", "unsafe")
            with zipfile.ZipFile(archive) as package, self.assertRaisesRegex(OSError, "unsafe path"):
                _safe_extract(package, root / "output")

    def test_parser_reports_invalid_gff(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "bad.gff3"
            path.write_text("chr\ttest\tCDS\t1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "expected 9"):
                list(parse_gff(path))

    def test_gff_parser_rejects_empty_core_fields_and_nonfinite_scores(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            empty_seqid = root / "empty-seqid.gff3"
            nonfinite_score = root / "nonfinite-score.gff3"
            empty_seqid.write_text(
                "\ttest\tCDS\t1\t3\t.\t+\t0\tID=x\n", encoding="utf-8"
            )
            nonfinite_score.write_text(
                "chr\ttest\tCDS\t1\t3\tnan\t+\t0\tID=x\n", encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "must not be empty"):
                list(parse_gff(empty_seqid))
            with self.assertRaisesRegex(ValueError, "score must be finite"):
                list(parse_gff(nonfinite_score))

    def test_fasta_parser_handles_wrapped_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "genome.fna"
            path.write_text(">one description\nACG\nTTA\n>two\nNNN\n", encoding="utf-8")
            self.assertEqual(parse_fasta(path), {"one": "ACGTTA", "two": "NNN"})

    def test_fasta_parser_rejects_empty_headers_and_records(self) -> None:
        cases = ((">\nACG\n", "empty FASTA identifier"), (">one\n", "contains no sequence"))
        for contents, message in cases:
            with self.subTest(contents=contents), tempfile.TemporaryDirectory() as temporary_directory:
                path = Path(temporary_directory) / "genome.fna"
                path.write_text(contents, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    parse_fasta(path)

    def test_cli_supports_direct_file_execution(self) -> None:
        cli_path = Path(__file__).parents[1] / "annostat" / "cli.py"
        result = subprocess.run(
            [sys.executable, str(cli_path), "--help"],
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Analyze bacterial GFF3", result.stdout)

    def test_cli_exposes_comparison_and_ncbi_help(self) -> None:
        cli_path = Path(__file__).parents[1] / "annostat" / "cli.py"
        comparison = subprocess.run(
            [sys.executable, str(cli_path), "compare", "--help"],
            check=False, capture_output=True, text=True,
        )
        fetch = subprocess.run(
            [sys.executable, str(cli_path), "fetch", "--help"],
            check=False, capture_output=True, text=True,
        )

        self.assertEqual(comparison.returncode, 0, comparison.stderr)
        self.assertIn("--genome LABEL FASTA GFF3", comparison.stdout)
        self.assertIn("--reference GCF_OR_GCA", comparison.stdout)
        self.assertEqual(fetch.returncode, 0, fetch.stderr)
        self.assertIn("official NCBI Datasets CLI", fetch.stdout)

    def test_cli_reports_version(self) -> None:
        cli_path = Path(__file__).parents[1] / "annostat" / "cli.py"
        result = subprocess.run(
            [sys.executable, str(cli_path), "--version"],
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertRegex(result.stdout, r"^Annostat \d+\.\d+\.\d+")

    def test_cli_rejects_an_inverted_cds_length_range(self) -> None:
        cli_path = Path(__file__).parents[1] / "annostat" / "cli.py"
        result = subprocess.run(
            [
                sys.executable,
                str(cli_path),
                "-f",
                "missing.fna",
                "-g",
                "missing.gff3",
                "--min-cds-length",
                "500",
                "--max-cds-length",
                "100",
            ],
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 2)
        self.assertIn("cannot be greater", result.stderr)

    def test_cli_prints_progress_summary_and_supports_quiet_mode(self) -> None:
        cli_path = Path(__file__).parents[1] / "annostat" / "cli.py"
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fasta = root / "genome.fna"
            gff = root / "genes.gff3"
            fasta.write_text(">chr\nATGAAATAA\n", encoding="utf-8")
            gff.write_text(
                "chr\ttest\tCDS\t1\t9\t.\t+\t0\tID=cds1;product=test;Dbxref=COG:J\n",
                encoding="utf-8",
            )
            command = [
                sys.executable, str(cli_path), "-f", str(fasta), "-g", str(gff),
                "-o", str(root / "results"),
            ]

            result = subprocess.run(command, check=False, capture_output=True, text=True)
            quiet_result = subprocess.run(
                command[:-1] + [str(root / "quiet-results"), "--quiet"],
                check=False,
                capture_output=True,
                text=True,
            )
            profile_result = subprocess.run(
                command[:-1] + [str(root / "profile-results"), "--profile"],
                check=False,
                capture_output=True,
                text=True,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("[1/5] Reading GFF3", result.stdout)
        self.assertIn("Analysis summary", result.stdout)
        self.assertIn("ATG/GTG/TTG starts", result.stdout)
        self.assertEqual(quiet_result.returncode, 0, quiet_result.stderr)
        self.assertEqual(quiet_result.stdout, "")
        self.assertEqual(profile_result.returncode, 0, profile_result.stderr)
        self.assertIn("Performance profile", profile_result.stdout)
        self.assertIn("Peak Python memory", profile_result.stdout)


if __name__ == "__main__":
    unittest.main()
