"""Focused regressions for the repaired bio-codon-usage surfaces."""

from __future__ import annotations

import contextlib
import io
import runpy
import sys
import unittest
import warnings
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqUtils import CodonAdaptationIndex


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from codon_utils import (  # noqa: E402
    CAIScoringError,
    CDSValidationError,
    calculate_cai_guarded,
    count_codons,
    optimize_dna_table_aware,
    validate_cds,
)
from rscu_analysis import calculate_rscu  # noqa: E402


class ValidationTests(unittest.TestCase):
    def test_strict_validates_full_cds(self):
        result = validate_cds("atggcttaa", table_id=1)
        self.assertEqual(result.sequence, "ATGGCTTAA")
        self.assertTrue(result.starts_with_allowed_codon)
        self.assertTrue(result.has_terminal_stop)
        self.assertEqual(result.discarded, ())

    def test_strict_rejects_empty_partial_shifted_ambiguous_and_internal_stop(self):
        cases = ("", "ATGAA", "AATGCTTAA", "ATGNNNTAA", "ATGTAAGCTTAA")
        for sequence in cases:
            with self.subTest(sequence=sequence), self.assertRaises(CDSValidationError):
                validate_cds(sequence, table_id=1)

    def test_permissive_discards_and_reports_without_shifting_triplets(self):
        result = validate_cds("ATGNNNGCTTAAAA", table_id=1, policy="permissive")
        self.assertEqual(result.sequence, "ATGGCTTAA")
        self.assertEqual(
            [(item.offset, item.text, item.reason) for item in result.discarded],
            [
                (3, "NNN", "codon contains non-ACGT symbols"),
                (12, "AA", "incomplete trailing codon"),
            ],
        )
        self.assertIn("offset 3", result.discard_summary())
        self.assertIn("offset 12", result.discard_summary())

    def test_table_two_tga_tgg_are_sense_and_aga_agg_are_stops(self):
        for sense, stop in (("TGA", "AGA"), ("TGG", "AGG"), ("TGA", "AGG"), ("TGG", "AGA")):
            with self.subTest(sense=sense, stop=stop):
                result = validate_cds(f"ATG{sense}{stop}", table_id=2)
                self.assertEqual(result.codons, ("ATG", sense, stop))

    def test_counts_and_rscu_share_validation(self):
        counts, count_check = count_codons("ATGCTTCTGTAA", table_id=1)
        rscu, rscu_check = calculate_rscu("ATGCTTCTGTAA", table_id=1)
        self.assertEqual(count_check, rscu_check)
        self.assertEqual(counts["CTT"], 1)
        self.assertAlmostEqual(rscu["CTT"] + rscu["CTG"], 6.0)


class CAITests(unittest.TestCase):
    def setUp(self):
        self.standard = CodonAdaptationIndex(
            [Seq("ATGGCTGCTGCTTAA"), Seq("ATGGCTGCTGCTTAA")],
            table=CodonTable.unambiguous_dna_by_id[1],
        )

    def test_guarded_scoring_reports_exclusions(self):
        result = calculate_cai_guarded(self.standard, "ATGGCTTAA", table_id=1)
        self.assertGreater(result.score, 0)
        self.assertEqual(result.included_codons, 2)
        self.assertEqual(result.excluded_codons[0][1], "ATG")

    def test_guard_prevents_zero_denominator(self):
        with self.assertRaisesRegex(CAIScoringError, "no included codons"):
            calculate_cai_guarded(
                self.standard,
                "ATGTGG",
                table_id=1,
                require_terminal_stop=False,
            )

    def test_biopython_185_semantics_are_pinned(self):
        stop_index = CodonAdaptationIndex([Seq("ATGGCTTAAGCTTAA")])
        self.assertLess(stop_index.calculate(Seq("GCTTAG")), stop_index.calculate(Seq("GCT")))
        self.assertAlmostEqual(CodonAdaptationIndex([Seq("GCT" * 10)])["GCC"], 0.05)
        tie_index = CodonAdaptationIndex([Seq("GCTGCCGCAGCG")])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            tie_index.optimize(Seq("GCT"), strict=False)
        self.assertEqual(caught, [])

    def test_table_two_optimization_preserves_all_target_codons(self):
        table2 = CodonTable.unambiguous_dna_by_id[2]
        cai = CodonAdaptationIndex([Seq("ATGTGATGAAGA")], table=table2)
        for sense, stop in (("TGA", "AGA"), ("TGG", "AGG"), ("TGA", "AGG"), ("TGG", "AGA")):
            query = Seq(f"ATG{sense}{stop}")
            with self.subTest(query=str(query)):
                optimized, source_check, optimized_check = optimize_dna_table_aware(
                    cai, query, table_id=2, strict=False
                )
                self.assertEqual(source_check.codons, ("ATG", sense, stop))
                self.assertEqual(
                    optimized.translate(table=2),
                    query.translate(table=2),
                )
                self.assertIn(optimized_check.codons[-1], {"AGA", "AGG"})

    def test_table_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "not constructed"):
            calculate_cai_guarded(self.standard, "ATGTGAAGA", table_id=2)


class SurfaceTests(unittest.TestCase):
    def test_scripts_are_import_safe(self):
        for script in ("basic_analysis.py", "rscu_analysis.py", "cai_optimization.py"):
            with self.subTest(script=script):
                stdout = io.StringIO()
                stderr = io.StringIO()
                with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                    runpy.run_path(str(SCRIPTS / script), run_name="not_main")
                self.assertEqual(stdout.getvalue(), "")
                self.assertEqual(stderr.getvalue(), "")

    def test_nonstandard_nc_helper_is_not_shipped_or_advertised(self):
        metrics = (ROOT / "references" / "metrics-and-methods.md").read_text(encoding="utf-8")
        skill = (ROOT / "SKILL.md").read_text(encoding="utf-8")
        self.assertNotIn("def effective_nc", metrics)
        self.assertNotIn("RSCU, and Nc", skill)
        self.assertIn("does not calculate Nc", skill)


if __name__ == "__main__":
    unittest.main(verbosity=2)
