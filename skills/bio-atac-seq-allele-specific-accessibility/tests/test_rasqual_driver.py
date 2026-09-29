import json
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_rasqual_features.py"


class RasqualDriverTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        samples = [f"sample{i}" for i in range(10)]
        with gzip.open(self.root / "cohort.vcf.gz", "wt", encoding="utf-8") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]) + "\n")

    def tearDown(self):
        self.temp.cleanup()

    def run_driver(self, feature_text, matrix_rows=3):
        matrix = bytes(10 * matrix_rows * 8)
        (self.root / "counts.bin").write_bytes(matrix)
        (self.root / "offsets.bin").write_bytes(matrix)
        features = self.root / "features.tsv"
        features.write_text(feature_text, encoding="utf-8")
        return subprocess.run(
            [
                sys.executable, str(SCRIPT), "--features", str(features),
                "--counts", str(self.root / "counts.bin"), "--offsets",
                str(self.root / "offsets.bin"), "--vcf",
                str(self.root / "cohort.vcf.gz"), "--samples", "10",
                "--matrix-rows", str(matrix_rows), "--output-dir",
                str(self.root / "out"), "--dry-run",
            ],
            text=True, capture_output=True, check=False,
        )

    def test_feature_rows_advance_and_zero_fsnp_row_keeps_matrix_index(self):
        completed = self.run_driver(
            "featureA\tchr1\t100\t200\t5\t2\n"
            "featureSkip\tchr1\t300\t400\t4\t0\n"
            "featureC\tchr2\t500\t600\t7\t3\n"
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        plan = json.loads(completed.stdout)
        self.assertEqual([item["index"] for item in plan["features"]], [1, 3])
        self.assertEqual(
            [item["rasqual"][item["rasqual"].index("-j") + 1] for item in plan["features"]],
            ["1", "3"],
        )
        self.assertEqual(plan["skipped_zero_feature_snps"], ["featureSkip"])
        self.assertIn("BH", plan["multiplicity"])

    def test_feature_row_cannot_exceed_declared_matrix(self):
        completed = self.run_driver(
            "featureA\tchr1\t100\t200\t5\t2\n"
            "featureB\tchr1\t300\t400\t4\t1\n",
            matrix_rows=1,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("exceeds declared matrix rows", completed.stderr)

    def test_feature_contract_rejects_inconsistent_counts(self):
        completed = self.run_driver("bad\tchr1\t100\t200\t2\t3\n")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("inconsistent SNP counts", completed.stderr)


if __name__ == "__main__":
    unittest.main()
