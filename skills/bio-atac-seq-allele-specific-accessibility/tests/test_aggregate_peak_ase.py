import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "aggregate_peak_ase.py"
ASE_HEADER = [
    "contig", "position", "variantID", "refAllele", "altAllele",
    "refCount", "altCount", "totalCount", "lowMAPQDepth", "lowBaseQDepth",
    "rawDepth", "otherBases", "improperPairs",
]


class AggregatePeakAseTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.peaks = self.root / "peaks.bed"
        self.peaks.write_text("chrAudit\t99\t250\tpeakA\n", encoding="utf-8")

    def tearDown(self):
        self.temp.cleanup()

    def write_counts(self, rows):
        path = self.root / "counts.tsv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(ASE_HEADER)
            writer.writerows(rows)
        return path

    def write_map(self, rows):
        path = self.root / "haplotypes.tsv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "contig", "position", "refAllele", "altAllele", "phaseSet",
                "haplotype1Allele",
            ])
            writer.writerows(rows)
        return path

    def run_helper(self, counts, haplotypes, peaks=None):
        output = self.root / "output.tsv"
        completed = subprocess.run(
            [
                sys.executable, str(SCRIPT), "--ase-counts", str(counts),
                "--peaks", str(peaks or self.peaks), "--haplotype-map",
                str(haplotypes), "--output", str(output),
            ],
            text=True, capture_output=True, check=False,
        )
        return completed, output

    def test_opposite_reference_labels_preserve_haplotype_direction(self):
        counts = self.write_counts([
            ["chrAudit", 100, "a", "A", "G", 90, 10, 100, 0, 0, 100, 0, 0],
            ["chrAudit", 200, "b", "C", "T", 10, 90, 100, 0, 0, 100, 0, 0],
        ])
        haplotypes = self.write_map([
            ["chrAudit", 100, "A", "G", "chrAudit:10", "REF"],
            ["chrAudit", 200, "C", "T", "chrAudit:10", "ALT"],
        ])
        completed, output = self.run_helper(counts, haplotypes)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        with output.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["haplotype1_count"], "180")
        self.assertEqual(rows[0]["haplotype2_count"], "20")
        self.assertEqual(rows[0]["status"], "significant_imbalance")

    def test_separate_phase_sets_are_not_pooled(self):
        counts = self.write_counts([
            ["chrAudit", 100, "a", "A", "G", 90, 10, 100, 0, 0, 100, 0, 0],
            ["chrAudit", 200, "b", "C", "T", 90, 10, 100, 0, 0, 100, 0, 0],
        ])
        haplotypes = self.write_map([
            ["chrAudit", 100, "A", "G", "chrAudit:10", "REF"],
            ["chrAudit", 200, "C", "T", "chrAudit:20", "REF"],
        ])
        completed, output = self.run_helper(counts, haplotypes)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        with output.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["status"] for row in rows}, {"underpowered_snp_count"})

    def test_disjoint_peaks_with_duplicate_display_names_are_not_pooled(self):
        counts = self.write_counts([
            ["chrAudit", 100, "a", "A", "G", 90, 10, 100, 0, 0, 100, 0, 0],
            ["chrAudit", 300, "b", "C", "T", 10, 90, 100, 0, 0, 100, 0, 0],
        ])
        haplotypes = self.write_map([
            ["chrAudit", 100, "A", "G", "chrAudit:10", "REF"],
            ["chrAudit", 300, "C", "T", "chrAudit:10", "ALT"],
        ])
        duplicate_names = self.root / "duplicate_names.bed"
        duplicate_names.write_text(
            "chrAudit\t99\t101\tdup\nchrAudit\t299\t301\tdup\n",
            encoding="utf-8",
        )

        completed, output = self.run_helper(counts, haplotypes, duplicate_names)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        with output.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        self.assertEqual([row["peak"] for row in rows], ["dup", "dup"])
        self.assertEqual([row["snp_count"] for row in rows], ["1", "1"])
        self.assertEqual(
            {row["status"] for row in rows},
            {"underpowered_snp_count"},
        )

    def test_valid_empty_cases_write_schema_correct_empty_table(self):
        counts = self.write_counts([
            ["chrAudit", 100, "a", "A", "G", 10, 10, 20, 0, 0, 20, 0, 0],
        ])
        haplotypes = self.write_map([
            ["chrAudit", 100, "A", "G", "chrAudit:10", "REF"],
        ])
        completed, output = self.run_helper(counts, haplotypes)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        with output.open(encoding="utf-8") as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 0)
        self.assertIn("all sites are below", completed.stdout)

        no_overlap = self.root / "no_overlap.bed"
        no_overlap.write_text("chrAudit\t1000\t1100\n", encoding="utf-8")
        counts = self.write_counts([
            ["chrAudit", 100, "a", "A", "G", 20, 10, 30, 0, 0, 30, 0, 0],
        ])
        completed, output = self.run_helper(counts, haplotypes, no_overlap)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        with output.open(encoding="utf-8") as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 0)
        self.assertIn("no depth-eligible sites overlap", completed.stdout)

    def test_malformed_schema_has_stable_actionable_error(self):
        counts = self.root / "bad.tsv"
        counts.write_text("contig\tposition\trefCount\nchrAudit\t100\t30\n", encoding="utf-8")
        haplotypes = self.write_map([
            ["chrAudit", 100, "A", "G", "chrAudit:10", "REF"],
        ])
        completed, output = self.run_helper(counts, haplotypes)
        self.assertEqual(completed.returncode, 2)
        self.assertIn("ASE table is missing required columns", completed.stderr)
        self.assertFalse(output.exists())

    def test_malformed_peak_schema_has_stable_actionable_error(self):
        counts = self.write_counts([
            ["chrAudit", 100, "a", "A", "G", 20, 10, 30, 0, 0, 30, 0, 0],
        ])
        haplotypes = self.write_map([
            ["chrAudit", 100, "A", "G", "chrAudit:10", "REF"],
        ])
        bad_peaks = self.root / "bad_peaks.bed"
        bad_peaks.write_text("chrAudit\tnot-an-integer\t200\n", encoding="utf-8")
        completed, output = self.run_helper(counts, haplotypes, bad_peaks)
        self.assertEqual(completed.returncode, 2)
        self.assertIn("peak start and end columns must contain integers", completed.stderr)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
