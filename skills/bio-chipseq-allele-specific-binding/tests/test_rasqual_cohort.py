import importlib.util
from pathlib import Path
import tempfile
import unittest


MODULE = Path(__file__).resolve().parents[1] / "scripts" / "rasqual_cohort.py"
SPEC = importlib.util.spec_from_file_location("rasqual_cohort", MODULE)
rasqual = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(rasqual)


class RasqualContractTests(unittest.TestCase):
    def test_bh_is_monotone_in_p_value_order(self):
        p_values = [0.04, 0.001, 0.03, 0.2]
        q_values = rasqual.benjamini_hochberg(p_values)
        ordered = sorted(zip(p_values, q_values))
        self.assertEqual(q_values, [0.05333333333333334, 0.004, 0.05333333333333334, 0.2])
        self.assertEqual([q for _, q in ordered], sorted(q for _, q in ordered))

    def test_parse_requires_25_columns_and_convergence(self):
        row = ["feature", *[str(i) for i in range(1, 25)]]
        row[13] = "0.5"  # column 14, phi
        row[22] = "0"    # column 23, convergence
        self.assertEqual(len(rasqual.parse_output("\t".join(row), "feature")), 1)
        row[22] = "1"
        with self.assertRaisesRegex(rasqual.ContractError, "convergence"):
            rasqual.parse_output("\t".join(row), "feature")

    def test_manifest_must_cover_each_feature_once(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.tsv"
            path.write_text(
                "\t".join(rasqual.MANIFEST_COLUMNS) + "\n" +
                "f1\t1:1-2\t1\t2\t1\t1\t1\t2\ttrue\n" +
                "f2\t1:3-4\t1\t2\t1\t1\t3\t4\ttrue\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(rasqual.ContractError, "exactly once"):
                rasqual.read_manifest(path, 2)


if __name__ == "__main__":
    unittest.main()
