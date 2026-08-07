from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from annostat.analysis import analyze_codons, analyze_features, cog_categories
from annostat.cli import run_analysis
from annostat.models import Feature
from annostat.parsers import parse_attributes, parse_fasta, parse_gff
from annostat.sequences import extract_cds_sequences, reverse_complement, translate_dna


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

            summary = run_analysis(fasta, gff, output, "tsv")

            self.assertEqual(summary["cds_count"], 1)
            expected = {
                "features.tsv", "cds_nucleotide.fasta", "cds_protein.fasta",
                "codon_usage.csv", "start_codons.csv", "cog_categories.csv",
                "summary.json", "cog_categories.svg", "cds_lengths.svg",
            }
            self.assertEqual({path.name for path in output.iterdir()}, expected)
            with (output / "features.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["ID"], "cds1")
            self.assertEqual(json.loads((output / "summary.json").read_text())["cds_count"], 1)

    def test_parser_reports_invalid_gff(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "bad.gff3"
            path.write_text("chr\ttest\tCDS\t1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "expected 9"):
                list(parse_gff(path))

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


if __name__ == "__main__":
    unittest.main()
