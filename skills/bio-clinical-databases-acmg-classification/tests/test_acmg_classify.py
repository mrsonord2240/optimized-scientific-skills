import importlib.util
import pathlib
import subprocess
import sys
import unittest
from unittest.mock import Mock, patch


SCRIPT = pathlib.Path(__file__).parents[1] / "scripts" / "acmg_classify.py"
FRAMEWORK = pathlib.Path(__file__).parents[1] / "references" / "acmg-framework.md"
SPEC = importlib.util.spec_from_file_location("acmg_classify", SCRIPT)
acmg = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(acmg)


class OddsPathTests(unittest.TestCase):
    def test_exact_brnich_boundaries(self):
        cases = {
            350.000001: "PS3_VeryStrong", 350: "PS3", 18.700001: "PS3",
            18.7: "PS3_Moderate", 4.300001: "PS3_Moderate",
            4.3: "PS3_Supporting", 2.100001: "PS3_Supporting", 2.1: None,
            0.48: None, 0.479999: "BS3_Supporting", 0.23: "BS3_Supporting",
            0.229999: "BS3_Moderate", 0.053: "BS3_Moderate", 0.052999: "BS3",
        }
        for value, expected in cases.items():
            with self.subTest(value=value):
                self.assertEqual(acmg.ps3_oddspath(value), expected)

    def test_nonpositive_rejected(self):
        for value in (0, -1):
            with self.assertRaises(ValueError):
                acmg.ps3_oddspath(value)

    def test_framework_prose_matches_brnich_boundaries(self):
        framework = FRAMEWORK.read_text(encoding="utf-8")
        self.assertIn(
            "OddsPath > 4.3 for Moderate and > 18.7 for Strong", framework
        )
        self.assertNotIn("OddsPath > 4.3 for Strong", framework)


class PVS1Tests(unittest.TestCase):
    def args(self, **overrides):
        values = dict(variant_type="nonsense", is_nmd_predicted=True,
                      coding_pct_removed=0.01, in_critical_region=False,
                      is_disease_relevant_transcript=True,
                      lof_mechanism_established=True)
        values.update(overrides)
        return values

    def test_prerequisites_gate_pvs1(self):
        self.assertIsNone(acmg.pvs1_decision_tree(**self.args(lof_mechanism_established=False)))
        self.assertIsNone(acmg.pvs1_decision_tree(**self.args(is_disease_relevant_transcript=False)))

    def test_nmd_and_ten_percent_boundary(self):
        self.assertEqual(acmg.pvs1_decision_tree(**self.args()), "PVS1_VeryStrong")
        self.assertEqual(acmg.pvs1_decision_tree(**self.args(is_nmd_predicted=False, coding_pct_removed=0.10)), "PVS1_Moderate")
        self.assertEqual(acmg.pvs1_decision_tree(**self.args(is_nmd_predicted=False, coding_pct_removed=0.100001)), "PVS1_Strong")

    def test_splice_requires_consequence_and_rescue_review(self):
        base = self.args(variant_type="splice_donor", is_nmd_predicted=True,
                         coding_pct_removed=0)
        self.assertEqual(acmg.pvs1_decision_tree(**base), acmg.PVS1_REVIEW_REQUIRED)
        self.assertEqual(acmg.pvs1_decision_tree(**base, splice_consequence="nmd"), acmg.PVS1_REVIEW_REQUIRED)
        self.assertEqual(acmg.pvs1_decision_tree(**base, splice_consequence="nmd", rescue_transcript_excluded=True), "PVS1_VeryStrong")
        self.assertIsNone(acmg.pvs1_decision_tree(
            **self.args(variant_type="splice_donor", is_nmd_predicted=False,
                        coding_pct_removed=0),
            splice_consequence="no_impact", rescue_transcript_excluded=True,
        ))

    def test_splice_rejects_contradictory_nmd_representations(self):
        for predicted, consequence in (
            (False, "nmd"), (True, "in_frame_disruptive"), (True, "no_impact")
        ):
            with self.subTest(predicted=predicted, consequence=consequence):
                with self.assertRaisesRegex(ValueError, "conflicts with splice_consequence"):
                    acmg.pvs1_decision_tree(
                        **self.args(variant_type="splice_acceptor",
                                    is_nmd_predicted=predicted, coding_pct_removed=0),
                        splice_consequence=consequence, rescue_transcript_excluded=True,
                    )

    def test_initiation_and_deletion_paths_remain_bounded(self):
        self.assertEqual(
            acmg.pvs1_decision_tree(**self.args(variant_type="initiation")),
            "PVS1_Moderate",
        )
        self.assertEqual(
            acmg.pvs1_decision_tree(**self.args(
                variant_type="multi_exon_del", is_nmd_predicted=False,
                coding_pct_removed=0.05, in_critical_region=False,
            )),
            "PVS1_Strong",
        )


class AlphaMissenseTests(unittest.TestCase):
    def test_bergquist_intervals(self):
        cases = {
            0.070: "BP4_3pt", 0.071: "BP4_Moderate", 0.099: "BP4_Moderate",
            0.100: "BP4_Supporting", 0.169: "BP4_Supporting", 0.170: None,
            0.791: None, 0.792: "PP3_Supporting", 0.905: "PP3_Supporting",
            0.906: "PP3_Moderate", 0.971: "PP3_Moderate", 0.972: "PP3_3pt",
            0.989: "PP3_3pt", 0.990: "PP3_Strong",
        }
        for score, expected in cases.items():
            with self.subTest(score=score):
                self.assertEqual(acmg.alphamissense_pp3_bp4(score), expected)

    def test_probability_domain(self):
        for score in (-0.001, 1.001):
            with self.assertRaises(ValueError):
                acmg.alphamissense_pp3_bp4(score)


class PejaverPredictorTests(unittest.TestCase):
    def test_generic_calibrator_retains_half_open_three_field_intervals(self):
        intervals = [(0.0, 0.5, "LOW"), (0.5, 1.0, "HIGH")]
        self.assertEqual(acmg.pejaver_calibrate(0.0, intervals), "LOW")
        self.assertEqual(acmg.pejaver_calibrate(0.5, intervals), "HIGH")
        self.assertIsNone(acmg.pejaver_calibrate(1.0, intervals))

    def test_revel_exact_endpoints_are_inclusive_on_the_benign_side(self):
        cases = {
            0.003: "BP4_VeryStrong", 0.016: "BP4_Strong",
            0.183: "BP4_Moderate", 0.290: "BP4_Supporting",
            0.644: "PP3_Supporting", 0.773: "PP3_Moderate",
            0.932: "PP3_Strong",
        }
        for score, expected in cases.items():
            with self.subTest(score=score):
                self.assertEqual(acmg.revel_pp3_bp4(score), expected)

    def test_bayesdel_exact_endpoints_are_inclusive_on_the_benign_side(self):
        cases = {
            -0.36: "BP4_Moderate", -0.18: "BP4_Supporting",
            0.13: "PP3_Supporting", 0.27: "PP3_Moderate",
            0.50: "PP3_Strong",
        }
        for score, expected in cases.items():
            with self.subTest(score=score):
                self.assertEqual(acmg.bayesdel_pp3_bp4(score), expected)


class ValidationTests(unittest.TestCase):
    def test_unknown_duplicate_and_wrong_type_rejected(self):
        with self.assertRaises(ValueError):
            acmg.tavtigian_classify(["NOT_A_CRITERION"])
        with self.assertRaises(ValueError):
            acmg.tavtigian_classify(["PS3", "PS3"])
        with self.assertRaises(TypeError):
            acmg.tavtigian_classify("PS3")

    def test_subsumption_and_notice(self):
        result = acmg.classify_with_subsumption(["PVS1_VeryStrong", "PP3_Strong", "PM2_Supporting"])
        self.assertEqual(result["criteria"], ["PVS1_VeryStrong", "PM2_Supporting"])
        self.assertTrue(result["review_required"])
        self.assertIn("not for patient diagnosis", result["clinical_use"])

    def test_current_rule_criteria_reject_retired_and_same_family_codes(self):
        for criteria in (
            ["PP5"], ["BP6"],
            ["PVS1_VeryStrong", "PVS1_Strong"],
            ["PS3", "PS3_Moderate"],
            ["PM2", "PM2_Supporting"],
            ["PP3_Supporting", "PP3_Strong"],
        ):
            with self.subTest(criteria=criteria):
                with self.assertRaises(ValueError):
                    acmg.tavtigian_classify(criteria)

    def test_opposing_computational_evidence_requires_conflict_review(self):
        with self.assertRaisesRegex(ValueError, "opposing PP3 and BP4"):
            acmg.tavtigian_classify(["PP3_Strong", "BP4_Strong"])

    def test_numeric_domains(self):
        for call in (
            lambda: acmg.spliceai_walker2023(-0.1),
            lambda: acmg.revel_pp3_bp4(1.1),
            lambda: acmg.whiffin_max_credible_af(-0.1),
            lambda: acmg.whiffin_max_credible_af(0.1, penetrance=0),
            lambda: acmg.bs1_ba1(-0.1, 1e-6),
        ):
            with self.assertRaises((TypeError, ValueError)):
                call()

    def test_standalone_demo_is_coherent_and_guarded(self):
        completed = subprocess.run(
            [sys.executable, "-B", str(SCRIPT)], capture_output=True, text=True, check=True
        )
        self.assertIn("NON-DIAGNOSTIC TRAINING EXAMPLE", completed.stdout)
        self.assertIn("disease=unresolved", completed.stdout)
        self.assertIn("VCEP=not checked", completed.stdout)
        self.assertIn("assay validity=unreviewed", completed.stdout)
        self.assertIn("PS3_Moderate", completed.stdout)
        self.assertNotIn("BP4_Supporting", completed.stdout)


class InterfaceTests(unittest.TestCase):
    @patch.object(acmg.requests, "get")
    def test_genebe_uses_coordinate_contract(self, get):
        response = Mock()
        response.json.return_value = {"variants": []}
        get.return_value = response
        self.assertEqual(acmg.genebe_api("chr17", 28364411, "c", "a"), {"variants": []})
        get.assert_called_once_with(
            "https://api.genebe.net/cloud/api-public/v1/variant",
            params={"chr": "17", "pos": 28364411, "ref": "C", "alt": "A", "genome": "hg38"},
            timeout=30.0,
        )
        response.raise_for_status.assert_called_once()

    @patch.object(acmg.requests, "get")
    def test_cspec_uses_versioned_gene_route(self, get):
        response = Mock()
        response.json.return_value = {"data": []}
        get.return_value = response
        acmg.cspec_gene_versions("gatm")
        get.assert_called_once_with(
            "https://cspec.genome.network/cspec/Gene/id/GATM/SequenceVariantInterpretation/version",
            timeout=30.0,
        )

    def test_interface_validation_precedes_request(self):
        for call in (
            lambda: acmg.genebe_api("chr0", 1, "A", "C"),
            lambda: acmg.genebe_api("1", -1, "A", "C"),
            lambda: acmg.cspec_gene_versions("../bad"),
        ):
            with self.assertRaises(ValueError):
                call()


class SomaticTests(unittest.TestCase):
    def tier(self, level, significance="oncogenic", same=False):
        return acmg.cancer_amp_tier(
            level, same_tumor_type=same, clinical_significance=significance,
            knowledgebase="CIViC evidence record CIViC-1", access_date="2026-09-28",
        )["tier"]

    def test_all_tiers_and_tier_iii_iv_distinction(self):
        self.assertEqual(self.tier("regulatory_approved", same=True), "Tier I-A")
        self.assertEqual(self.tier("professional_guideline", same=True), "Tier I-B")
        self.assertEqual(self.tier("clinical_trial"), "Tier II-C")
        self.assertEqual(self.tier("preclinical"), "Tier II-D")
        self.assertEqual(self.tier("none", "uncertain"), "Tier III")
        self.assertEqual(self.tier("none", "benign"), "Tier IV")

    def test_invalid_or_unproven_somatic_input_stops(self):
        with self.assertRaises(ValueError):
            self.tier(-1)
        with self.assertRaises(ValueError):
            self.tier("none", "oncogenic")
        with self.assertRaises(ValueError):
            acmg.cancer_amp_tier("none", same_tumor_type=False,
                                 clinical_significance="uncertain", knowledgebase="",
                                 access_date="2026-09-28")
        for invalid_date in ("2026-99-99", "2025-02-29", "2026-9-28"):
            with self.subTest(access_date=invalid_date):
                with self.assertRaises(ValueError):
                    acmg.cancer_amp_tier(
                        "none", same_tumor_type=False,
                        clinical_significance="uncertain",
                        knowledgebase="CIViC evidence record CIViC-1",
                        access_date=invalid_date,
                    )


if __name__ == "__main__":
    unittest.main()
