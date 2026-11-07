#!/usr/bin/env python3
"""Focused regressions for the executable cfDNA workflow.

Set ``CFDNA_FIXTURE_ROOT`` to a directory containing ``raw.unmapped.bam``,
``reference.fa`` plus its bwa/FASTA/dictionary sidecars, and ``flag-qc.bam``.
The bounded tooling fixture prepared for this skill has exactly that layout.
"""

from __future__ import annotations

import ast
import importlib.util
import os
from pathlib import Path
import shutil
import tempfile
import unittest

import pysam


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "preprocess_cfdna.py"
SPEC = importlib.util.spec_from_file_location("preprocess_cfdna", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class ValidationTests(unittest.TestCase):
    def test_threads_and_umi_layout_are_validated(self):
        for value in (True, 0, -1, 257, 1.5):
            with self.subTest(value=value), self.assertRaises((TypeError, ValueError)):
                MODULE._positive_threads(value)
        self.assertEqual(MODULE._positive_threads(1), 1)

        with self.assertRaises(ValueError):
            MODULE._validated_umi_layout("6M11S+T", ("ZA",))
        with self.assertRaises(ValueError):
            MODULE._validated_umi_layout("6M11S+T", ("ZA", "ZA"))
        with self.assertRaises(ValueError):
            MODULE._validated_umi_layout("6M11S+T", ("too-long", "ZB"))

    def test_no_shell_true_call_exists(self):
        tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
        shell_true = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                shell_true.extend(
                    keyword
                    for keyword in node.keywords
                    if keyword.arg == "shell"
                    and isinstance(keyword.value, ast.Constant)
                    and keyword.value.value is True
                )
        self.assertEqual(shell_true, [])


@unittest.skipUnless(os.environ.get("CFDNA_FIXTURE_ROOT"), "set CFDNA_FIXTURE_ROOT")
class LiveWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = Path(os.environ["CFDNA_FIXTURE_ROOT"])
        cls.raw = cls.fixture / "raw.unmapped.bam"
        cls.reference = cls.fixture / "reference.fa"
        cls.flag_qc = cls.fixture / "flag-qc.bam"
        for path in (cls.raw, cls.reference, cls.flag_qc):
            if not path.is_file():
                raise RuntimeError(f"missing fixture: {path}")

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cfdna-regression-")
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def _run_and_assert(self, root: Path, *, duplex: bool = False) -> Path:
        output = root / "final.bam"
        returned = MODULE.preprocess_cfdna(
            self.raw,
            output,
            self.reference,
            duplex=duplex,
            threads=2,
        )
        self.assertEqual(returned, output)
        self.assertTrue(output.is_file())
        self.assertTrue(Path(f"{output}.bai").is_file())
        with pysam.AlignmentFile(output, "rb") as bam:
            self.assertEqual(bam.header.to_dict()["HD"]["SO"], "coordinate")
            self.assertGreater(sum(1 for _ in bam.fetch(until_eof=True)), 0)
        with pysam.AlignmentFile(root / "final_umis.bam", "rb", check_sq=False) as bam:
            records = list(bam.fetch(until_eof=True))
        self.assertGreater(len(records), 0)
        self.assertTrue(all(read.has_tag("RX") for read in records))
        self.assertTrue(all(read.has_tag("ZA") or read.has_tag("ZB") for read in records))
        return output

    def test_simplex_and_reciprocal_duplex(self):
        self._run_and_assert(self.root / "simplex", duplex=False)
        self._run_and_assert(self.root / "duplex", duplex=True)

    def test_space_and_metacharacter_paths_are_literal(self):
        unusual = self.root / "space ; dollar $ and [brackets]"
        unusual.mkdir()
        for source in self.fixture.glob("reference.fa*"):
            shutil.copy2(source, unusual / source.name)
        shutil.copy2(self.fixture / "reference.dict", unusual / "reference.dict")
        shutil.copy2(self.raw, unusual / "raw ; input $.bam")
        output = unusual / "output ; $ [safe]" / "final.bam"
        MODULE.preprocess_cfdna(
            unusual / "raw ; input $.bam",
            output,
            unusual / "reference.fa",
            threads=2,
        )
        with pysam.AlignmentFile(output, "rb") as bam:
            self.assertGreater(sum(1 for _ in bam.fetch(until_eof=True)), 0)
        self.assertFalse((unusual / "safe").exists())

    def test_fragment_qc_flags_bounds_and_empty_input(self):
        result = MODULE.insert_size_qc(self.flag_qc)
        counts = result["filtered_counts"]
        self.assertGreater(counts["secondary"], 0)
        self.assertGreater(counts["supplementary"], 0)
        self.assertGreater(counts["duplicate"], 0)
        self.assertGreater(counts["qc_fail"], 0)
        self.assertEqual(result["n"], counts["accepted"])
        self.assertEqual(
            result["interpretation_context"],
            [
                "library_chemistry",
                "collection_and_plasma_processing",
                "physical_or_in_silico_size_selection",
            ],
        )
        for invalid in (0, -1, True, 1.5):
            with self.subTest(invalid=invalid), self.assertRaises((TypeError, ValueError)):
                MODULE.insert_size_qc(self.flag_qc, max_size=invalid)

        empty = self.root / "empty.bam"
        with pysam.AlignmentFile(self.flag_qc, "rb") as source:
            with pysam.AlignmentFile(empty, "wb", header=source.header):
                pass
        empty_result = MODULE.insert_size_qc(empty)
        self.assertEqual(empty_result["n"], 0)
        self.assertEqual(empty_result["filtered_counts"]["input_records"], 0)

    def test_ten_consecutive_calls_are_stable(self):
        record_counts = []
        for index in range(10):
            output = self._run_and_assert(self.root / f"call-{index + 1:02d}")
            with pysam.AlignmentFile(output, "rb") as bam:
                record_counts.append(sum(1 for _ in bam.fetch(until_eof=True)))
        self.assertEqual(record_counts, [record_counts[0]] * 10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
