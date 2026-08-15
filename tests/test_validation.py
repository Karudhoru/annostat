"""Rigorous tests for validation, multipart CDS handling, and aggregation."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from annostat.analysis import analyze_features
from annostat.cohort import build_cohort, write_cohort
from annostat.parsers import parse_fasta, parse_gff
from annostat.qc import sequence_quality_findings
from annostat.sequences import (
    extract_cds_sequences,
    recognized_start_codons,
    resolve_genetic_code,
    reverse_complement,
    translate_dna,
)
from annostat.validation import RULES, canonical_json, validate_annotation, write_validation


class ValidationTests(unittest.TestCase):
    """Exercise deterministic rules with minimal, auditable fixtures."""

    def _write_pair(self, root: Path, fasta: str, gff: str) -> tuple[Path, Path]:
        fasta_path = root / "genome.fna"
        gff_path = root / "annotation.gff3"
        fasta_path.write_text(fasta, encoding="utf-8")
        gff_path.write_text(gff, encoding="utf-8")
        return fasta_path, gff_path

    def test_valid_pair_passes_with_hashes_and_rule_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fasta, gff = self._write_pair(
                Path(temporary_directory),
                ">chr\nATGAAATAA\n",
                "##gff-version 3\n"
                "chr\ttest\tgene\t1\t9\t.\t+\t.\tID=gene1\n"
                "chr\ttest\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=gene1;product=test\n",
            )

            result = validate_annotation(fasta, gff)

            self.assertTrue(result["valid"])
            self.assertEqual(result["severity_counts"], {"error": 0, "warning": 0, "info": 0})
            self.assertRegex(result["input_sha256"]["fasta"], r"^[0-9a-f]{64}$")
            self.assertRegex(result["scientific_fingerprint"], r"^[0-9a-f]{64}$")
            self.assertEqual({rule["rule_id"] for rule in result["rules"]}, set(RULES))
            self.assertTrue(all(rule.source_url.startswith("https://") for rule in RULES.values()))

    def test_cross_file_relationship_alphabet_and_phase_errors_are_separate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fasta, gff = self._write_pair(
                Path(temporary_directory),
                ">chr\nATGZAA\n",
                "chr\ttest\tCDS\t1\t12\t.\t+\t.\tID=bad;Parent=missing\n"
                "absent\ttest\tgene\t1\t3\t.\t+\t.\tID=other\n",
            )

            result = validate_annotation(fasta, gff)
            observed = {finding["rule_id"] for finding in result["findings"]}

            self.assertFalse(result["valid"])
            self.assertEqual(
                observed,
                {
                    "GFF_VERSION", "FASTA_ALPHABET", "COORDINATE_BOUNDS",
                    "SEQID_MATCH", "PARENT_EXISTS", "CDS_PHASE",
                },
            )
            self.assertEqual(result["severity_counts"]["error"], 5)
            self.assertEqual(result["severity_counts"]["warning"], 1)

    def test_multipart_plus_strand_is_counted_and_translated_once(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fasta, gff = self._write_pair(
                Path(temporary_directory),
                ">chr\nATGACCCCCAATAA\n",
                "##gff-version 3\n"
                "chr\ttest\tCDS\t1\t4\t.\t+\t0\tID=joined;product=test\n"
                "chr\ttest\tCDS\t10\t14\t.\t+\t2\tID=joined;product=test\n",
            )
            features = list(parse_gff(gff))

            result = validate_annotation(fasta, gff)
            records = extract_cds_sequences(features, parse_fasta(fasta))
            analysis = analyze_features(features)

            self.assertTrue(result["valid"])
            self.assertEqual(result["record_counts"]["cds_rows"], 2)
            self.assertEqual(analysis["gff_row_count"], 2)
            self.assertEqual(analysis["cds_count"], 1)
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].nucleotide, "ATGAAATAA")
            self.assertEqual(records[0].protein, "MK")
            self.assertEqual(len(records[0].segments), 2)

    def test_multipart_minus_strand_uses_biological_segment_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fasta, gff = self._write_pair(
                Path(temporary_directory),
                ">chr\nTTATTCCCCTCAT\n",
                "##gff-version 3\n"
                "chr\ttest\tCDS\t1\t5\t.\t-\t2\tID=joined\n"
                "chr\ttest\tCDS\t10\t13\t.\t-\t0\tID=joined\n",
            )

            records = extract_cds_sequences(list(parse_gff(gff)), parse_fasta(fasta))

            self.assertTrue(validate_annotation(fasta, gff)["valid"])
            self.assertEqual(records[0].nucleotide, "ATGAAATAA")
            self.assertEqual(records[0].protein, "MK")

    def test_inconsistent_multipart_phase_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fasta, gff = self._write_pair(
                Path(temporary_directory),
                ">chr\nATGACCCCCAATAA\n",
                "##gff-version 3\n"
                "chr\ttest\tCDS\t1\t4\t.\t+\t0\tID=joined\n"
                "chr\ttest\tCDS\t10\t14\t.\t+\t0\tID=joined\n",
            )

            result = validate_annotation(fasta, gff)

            self.assertFalse(result["valid"])
            phase = next(item for item in result["findings"] if item["rule_id"] == "CDS_PHASE")
            self.assertEqual(phase["observed"], "0")
            self.assertEqual(phase["expected"], "2")

    def test_circular_virtual_coordinates_require_a_real_start_coordinate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fasta, gff = self._write_pair(
                Path(temporary_directory),
                ">chr\nATGAAATAA\n",
                "##gff-version 3\n"
                "chr\ttest\tregion\t1\t9\t.\t+\t.\tID=region;Is_circular=true\n"
                "chr\ttest\tCDS\t10\t15\t.\t+\t0\tID=outside\n",
            )

            result = validate_annotation(fasta, gff)

            self.assertFalse(result["valid"])
            bounds = [
                item for item in result["findings"]
                if item["rule_id"] == "COORDINATE_BOUNDS"
            ]
            self.assertEqual(len(bounds), 1)
            self.assertIn("starts beyond", bounds[0]["message"])

    def test_reverse_complement_is_an_involution_for_iupac_dna(self) -> None:
        sequence = "ACGTRYSWKMBDHVN"
        self.assertEqual(reverse_complement(reverse_complement(sequence)), sequence)

    def test_supported_prokaryotic_genetic_codes_translate_tga_correctly(self) -> None:
        sequence = "ATGTGATAA"
        self.assertEqual(translate_dna(sequence, genetic_code=11), "M*")
        self.assertEqual(translate_dna(sequence, genetic_code=4), "MW")
        self.assertEqual(translate_dna(sequence, genetic_code=25), "MG")
        self.assertEqual(translate_dna("TTATAA", genetic_code=4), "M")
        self.assertEqual(translate_dna("TTATAA", genetic_code=11), "L")
        with self.assertRaisesRegex(ValueError, "unsupported NCBI genetic code"):
            translate_dna(sequence, genetic_code=1)

    def test_translation_table_is_validated_and_resolved_from_gff3(self) -> None:
        from annostat.models import Feature

        table_25 = Feature(
            "chr", "test", "CDS", 1, 9, None, "+", 0,
            {"ID": "code25", "transl_table": "25"},
        )
        self.assertEqual(resolve_genetic_code([table_25]), 25)
        self.assertEqual(resolve_genetic_code([], None), 11)
        self.assertEqual(resolve_genetic_code([], 4), 4)
        self.assertEqual(recognized_start_codons(11), {
            "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG",
        })
        with self.assertRaisesRegex(ValueError, "conflicts with GFF3"):
            resolve_genetic_code([table_25], 11)

        region_table = Feature(
            "chr", "test", "region", 1, 9, None, "+", None,
            {"ID": "region", "transl_table": "4"},
        )
        self.assertEqual(resolve_genetic_code([region_table]), 4)

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fasta, gff = self._write_pair(
                root,
                ">chr\nATGTGATAA\n",
                "##gff-version 3\n"
                "chr\ttest\tCDS\t1\t9\t.\t+\t0\tID=one;transl_table=25\n"
                "chr\ttest\tCDS\t1\t9\t.\t+\t0\tID=two;transl_table=11\n",
            )
            result = validate_annotation(fasta, gff)

        translation_findings = [
            finding for finding in result["findings"]
            if finding["rule_id"] == "TRANSLATION_TABLE"
        ]
        self.assertFalse(result["valid"])
        self.assertEqual(len(translation_findings), 1)
        self.assertEqual(translation_findings[0]["observed"], "11, 25")

        for raw_code, expected_message in (
            ("not-a-number", "not an integer"),
            ("1", "unsupported translation table"),
        ):
            with self.subTest(transl_table=raw_code), tempfile.TemporaryDirectory() as directory:
                fasta, gff = self._write_pair(
                    Path(directory),
                    ">chr\nATGTGATAA\n",
                    "##gff-version 3\n"
                    f"chr\ttest\tCDS\t1\t9\t.\t+\t0\tID=one;transl_table={raw_code}\n",
                )
                malformed_result = validate_annotation(fasta, gff)
                finding = next(
                    item for item in malformed_result["findings"]
                    if item["rule_id"] == "TRANSLATION_TABLE"
                )
                self.assertIn(expected_message, finding["message"])
                self.assertFalse(malformed_result["valid"])

    def test_linear_contig_boundary_suppresses_only_boundary_start_assumption(self) -> None:
        from annostat.models import Feature

        feature = Feature("chr", "test", "CDS", 1, 9, None, "+", 0, {"ID": "edge"})
        record = extract_cds_sequences([feature], {"chr": "TTTAAATAA" + "A" * 11})[0]

        linear = sequence_quality_findings(record, sequence_length=20, circular=False)
        circular = sequence_quality_findings(record, sequence_length=20, circular=True)

        self.assertNotIn("unrecognized_start_codon", {item.issue_type for item in linear})
        self.assertIn("unrecognized_start_codon", {item.issue_type for item in circular})

    def test_partial_boundaries_suppress_only_the_affected_cds_end(self) -> None:
        from annostat.models import Feature

        three_prime_partial = Feature(
            "chr", "test", "CDS", 1, 9, None, "+", 0,
            {"ID": "three", "partial": "true", "end_range": "9,."},
        )
        five_prime_partial = Feature(
            "chr", "test", "CDS", 1, 9, None, "+", 0,
            {"ID": "five", "partial": "true", "start_range": ".,1"},
        )
        bad_start = extract_cds_sequences(
            [three_prime_partial], {"chr": "CCCAAATAA"}
        )[0]
        bad_stop = extract_cds_sequences(
            [five_prime_partial], {"chr": "ATGAAAAAA"}
        )[0]
        partial_alternative_start = extract_cds_sequences(
            [five_prime_partial], {"chr": "TTGAAATAA"}
        )[0]

        self.assertIn(
            "unrecognized_start_codon",
            {
                item.issue_type for item in sequence_quality_findings(
                    bad_start, sequence_length=20, circular=True
                )
            },
        )
        self.assertIn(
            "missing_stop_codon",
            {
                item.issue_type for item in sequence_quality_findings(
                    bad_stop, sequence_length=20, circular=True
                )
            },
        )
        self.assertEqual(partial_alternative_start.protein, "LK")

    def test_position_specific_selenocysteine_exception_is_applied(self) -> None:
        from annostat.models import Feature

        feature = Feature(
            "chr", "test", "CDS", 1, 9, None, "+", 0,
            {"ID": "sec", "transl_except": "(pos:4..6,aa:Sec)"},
        )
        record = extract_cds_sequences([feature], {"chr": "ATGTGATAA"})[0]

        self.assertEqual(record.protein, "MU")
        self.assertEqual(record.translation_exception_indices, (1,))
        self.assertNotIn(
            "internal_stop_codon",
            {item.issue_type for item in sequence_quality_findings(record)},
        )

    def test_other_translation_exception_becomes_unknown_amino_acid(self) -> None:
        from annostat.models import Feature

        feature = Feature(
            "chr", "test", "CDS", 1, 9, None, "+", 0,
            {"ID": "other", "transl_except": "(pos:4..6,aa:OTHER)"},
        )
        record = extract_cds_sequences([feature], {"chr": "ATGTGATAA"})[0]

        self.assertEqual(record.protein, "MX")
        self.assertEqual(record.translation_exception_indices, (1,))
        self.assertNotIn(
            "internal_stop_codon",
            {item.issue_type for item in sequence_quality_findings(record)},
        )

    def test_terminal_translation_exception_is_not_lost_with_stop_trimming(self) -> None:
        from annostat.models import Feature

        feature = Feature(
            "chr", "test", "CDS", 1, 6, None, "+", 0,
            {"ID": "terminal-sec", "transl_except": "(pos:4..6,aa:Sec)"},
        )
        record = extract_cds_sequences([feature], {"chr": "ATGTGA"})[0]

        self.assertEqual(record.protein, "MU")
        self.assertEqual(record.translation_exception_indices, (1,))

    def test_validation_outputs_are_byte_reproducible(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fasta, gff = self._write_pair(
                root,
                ">chr\nATGAAATAA\n",
                "##gff-version 3\nchr\ttest\tCDS\t1\t9\t.\t+\t0\tID=cds1\n",
            )
            result = validate_annotation(fasta, gff)
            first = root / "first"
            second = root / "second"

            write_validation(first, result)
            write_validation(second, validate_annotation(fasta, gff))

            self.assertEqual(canonical_json(result), canonical_json(validate_annotation(fasta, gff)))
            self.assertEqual(
                (first / "validation.json").read_bytes(),
                (second / "validation.json").read_bytes(),
            )
            self.assertEqual(
                (first / "validation.tsv").read_bytes(),
                (second / "validation.tsv").read_bytes(),
            )

    def test_validate_cli_has_ci_friendly_exit_status(self) -> None:
        cli = Path(__file__).parents[1] / "annostat" / "cli.py"
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fasta, gff = self._write_pair(
                root,
                ">chr\nATGAAATAA\n",
                "##gff-version 3\nchr\ttest\tCDS\t1\t9\t.\t+\t0\tID=cds1\n",
            )
            passed = subprocess.run(
                [sys.executable, str(cli), "validate", "-f", str(fasta), "-g", str(gff),
                 "-o", str(root / "pass")],
                capture_output=True, text=True, check=False,
            )
            self.assertFalse((root / "pass").exists())
            self.assertIn("No validation issues found", passed.stdout)
            gff.write_text(
                "chr\ttest\tCDS\t1\t9\t.\t+\t0\tID=cds1\n",
                encoding="utf-8",
            )
            warned = subprocess.run(
                [sys.executable, str(cli), "validate", "-f", str(fasta), "-g", str(gff),
                 "-o", str(root / "warning")],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(warned.returncode, 0, warned.stderr)
            self.assertTrue((root / "warning" / "validation.json").is_file())
            self.assertTrue((root / "warning" / "validation.tsv").is_file())
            gff.write_text("bad\n", encoding="utf-8")
            failed = subprocess.run(
                [sys.executable, str(cli), "validate", "-f", str(fasta), "-g", str(gff),
                 "-o", str(root / "fail")],
                capture_output=True, text=True, check=False,
            )
            self.assertTrue((root / "fail" / "validation.json").is_file())
            self.assertTrue((root / "fail" / "validation.tsv").is_file())
            ignored = subprocess.run(
                [sys.executable, str(cli), "validate", "-f", str(fasta), "-g", str(gff),
                 "-o", str(root / "ignore"), "--fail-on", "never"],
                capture_output=True, text=True, check=False,
            )
            self.assertTrue((root / "ignore" / "validation.json").is_file())
            self.assertTrue((root / "ignore" / "validation.tsv").is_file())

        self.assertEqual(passed.returncode, 0, passed.stderr)
        self.assertIn("PASS", passed.stdout)
        self.assertEqual(failed.returncode, 1, failed.stderr)
        self.assertIn("FAIL", failed.stdout)
        self.assertEqual(ignored.returncode, 0, ignored.stderr)


class CohortTests(unittest.TestCase):
    """Verify normalized cohort reporting and missing-data semantics."""

    def _summary(self, sample_dir: Path, *, cog: bool, version: str = "1.0.0") -> None:
        sample_dir.mkdir()
        summary = {
            "annostat_version": version,
            "genome_length": 1000,
            "sequence_ids": ["chr"],
            "cds_count": 10,
            "rna_counts": {"tRNA": 2},
            "hypothetical_cds_count": 2,
            "cds_with_cog_count": 5 if cog else 0,
            "cog_data_available": cog,
            "quality_control": {
                "genome_gc_percent": 50.0,
                "coding_density_percent": 85.0,
                "severity_counts": {"warning": 1},
            },
            "validation": {
                "valid": True,
                "severity_counts": {"error": 0, "warning": 0},
                "input_sha256": {"fasta": "a" * 64, "gff3": "b" * 64},
            },
        }
        (sample_dir / "summary.json").write_text(json.dumps(summary), encoding="utf-8")

    def test_cohort_is_sorted_deterministic_and_preserves_missing_cog(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self._summary(root / "zeta", cog=False)
            self._summary(root / "alpha", cog=True)

            cohort = build_cohort([root])
            first = root / "cohort-one"
            second = root / "cohort-two"
            write_cohort(first, cohort)
            write_cohort(second, build_cohort([root / "zeta", root / "alpha"]))

            self.assertEqual([row["sample"] for row in cohort["samples"]], ["alpha", "zeta"])
            self.assertIsNone(cohort["samples"][1]["cog_coverage_percent"])
            for name in ("cohort.json", "cohort.tsv", "cohort.html"):
                self.assertEqual((first / name).read_bytes(), (second / name).read_bytes())
            self.assertIn("NA", (first / "cohort.tsv").read_text(encoding="utf-8"))
            self.assertNotIn("<script", (first / "cohort.html").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
