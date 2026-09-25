import importlib.util
import pathlib
import sys
import types
import unittest


ROOT = pathlib.Path(__file__).parents[3]
SCRIPT = ROOT / "skills" / "bio-vcf-statistics" / "examples" / "vcf_stats.py"


class Variant:
    def __init__(self, ref, alts, *, qual, is_snp=True, is_indel=False, filtered=False):
        self.REF = ref
        self.ALT = alts
        self.QUAL = qual
        self.is_snp = is_snp
        self.is_indel = is_indel
        self.FILTER = "LowQual" if filtered else None


class FakeVCF:
    variants = []

    def __init__(self, _path):
        self.closed = False

    def __iter__(self):
        return iter(self.variants)

    def close(self):
        self.closed = True


def load_example():
    fake_cyvcf2 = types.ModuleType("cyvcf2")
    fake_cyvcf2.VCF = FakeVCF
    previous = sys.modules.get("cyvcf2")
    sys.modules["cyvcf2"] = fake_cyvcf2
    try:
        spec = importlib.util.spec_from_file_location("vcf_stats_example", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        if previous is None:
            sys.modules.pop("cyvcf2", None)
        else:
            sys.modules["cyvcf2"] = previous


class VcfStatsRegressionTests(unittest.TestCase):
    def test_counts_every_multiallelic_snp_alt_for_titv(self):
        FakeVCF.variants = [Variant("A", ["G", "C"], qual=50.0)]

        stats = load_example().calculate_stats("fixture.vcf")

        self.assertEqual(stats["snps"], 1)
        self.assertEqual(stats["transitions"], 1)
        self.assertEqual(stats["transversions"], 1)

    def test_zero_quality_is_observed_but_missing_quality_is_not(self):
        FakeVCF.variants = [
            Variant("A", ["G"], qual=0.0),
            Variant("C", ["T"], qual=20.0),
            Variant("G", ["A"], qual=None),
        ]

        stats = load_example().calculate_stats("fixture.vcf")

        self.assertEqual(stats["qual_count"], 2)
        self.assertEqual(stats["qual_sum"], 20.0)


if __name__ == "__main__":
    unittest.main()
