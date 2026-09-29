#!/usr/bin/env python3
"""Focused regressions for BATCH-001 through BATCH-007."""

from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.dont_write_bytecode = True

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


SKILL_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_ROOT / "scripts"))
import batch_process as batch  # noqa: E402
import pyfastx_index as pyfastx_example  # noqa: E402


def write_fasta(path, records):
    return SeqIO.write(
        (SeqRecord(Seq(sequence), id=record_id, description="") for record_id, sequence in records),
        path,
        "fasta",
    )


class BatchProcessRegressions(unittest.TestCase):
    def test_prefix_split_rejects_escape_and_preserves_sibling(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "hostile.fasta"
            output = root / "out"
            sibling = root / "escaped.fasta"
            sibling.write_text("sentinel\n", encoding="utf-8")
            write_fasta(source, [("../escaped_record", "ACGT")])

            with self.assertRaisesRegex(ValueError, "unsafe record-id prefix"):
                batch.split_by_prefix(source, "fasta", output, max_open_files=2)

            self.assertEqual(sibling.read_text(encoding="utf-8"), "sentinel\n")
            manifest = json.loads((output / "split-prefix-manifest.json").read_text())
            self.assertEqual(manifest["status"], "partial")
            self.assertEqual(manifest["records_written"], 0)
            self.assertEqual(manifest["error"]["type"], "ValueError")

    def test_prefix_split_never_overwrites_existing_target(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.fasta"
            output = root / "out"
            output.mkdir()
            target = output / "sample.fasta"
            target.write_text("sentinel\n", encoding="utf-8")
            write_fasta(source, [("sample_1", "ACGT")])

            with self.assertRaises(FileExistsError):
                batch.split_by_prefix(source, "fasta", output, max_open_files=2)
            self.assertEqual(target.read_text(encoding="utf-8"), "sentinel\n")

    def test_prefix_split_bounds_handles_and_writes_every_record(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "many.fasta"
            output = root / "out"
            expected = 2048
            write_fasta(source, ((f"p{i:04d}_1", "ACGT") for i in range(expected)))
            before = len(os.listdir("/proc/self/fd")) if Path("/proc/self/fd").is_dir() else None

            manifest = batch.split_by_prefix(source, "fasta", output, max_open_files=8)

            after = len(os.listdir("/proc/self/fd")) if before is not None else None
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["records_written"], expected)
            self.assertEqual(len(manifest["files"]), expected)
            self.assertEqual(sum(item["records"] for item in manifest["files"]), expected)
            if before is not None:
                self.assertLessEqual(after, before + 1)

    def test_prefix_split_reopens_evicted_group_without_truncating(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "interleaved.fasta"
            output = root / "out"
            write_fasta(
                source,
                [("a_1", "A"), ("b_1", "CC"), ("c_1", "GGG"), ("a_2", "TTTT")],
            )
            manifest = batch.split_by_prefix(source, "fasta", output, max_open_files=2)
            counts = {item["prefix"]: item["records"] for item in manifest["files"]}
            self.assertEqual(counts, {"a": 2, "b": 1, "c": 1})
            with (output / "a.fasta").open("r", encoding="utf-8") as handle:
                self.assertEqual([record.id for record in SeqIO.parse(handle, "fasta")], ["a_1", "a_2"])

    def test_prefix_split_rejects_case_insensitive_collision(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "collision.fasta"
            output = root / "out"
            write_fasta(source, [("Sample_1", "AC"), ("sample_2", "GT")])
            with self.assertRaisesRegex(ValueError, "case-insensitive filesystems"):
                batch.split_by_prefix(source, "fasta", output, max_open_files=2)
            manifest = json.loads((output / "split-prefix-manifest.json").read_text())
            self.assertEqual(manifest["status"], "partial")
            self.assertEqual(manifest["records_written"], 1)

    def test_invalid_chunk_sizes_fail_before_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.fasta"
            write_fasta(source, [("one", "ACGT")])
            for invalid in (0, -1, 1.5, True):
                with self.subTest(invalid=invalid):
                    with self.assertRaisesRegex(ValueError, "records_per_file must be a positive integer"):
                        batch.split_by_count(source, "fasta", invalid, root / "chunk")
            self.assertEqual(list(root.glob("chunk_*")), [])
            outputs = batch.split_by_count(source, "fasta", 1, root / "chunk")
            self.assertEqual(len(outputs), 1)
            with Path(outputs[0]).open("r", encoding="utf-8") as handle:
                self.assertEqual(sum(1 for _ in SeqIO.parse(handle, "fasta")), 1)

    def test_empty_summary_is_header_only_with_explicit_schema(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "summary.csv"
            self.assertEqual(batch.summarize_files([], "fasta", output), [])
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            self.assertEqual(rows, [list(batch.SUMMARY_FIELDS)])

    def test_stable_discovery_and_summary_ignore_creation_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            left, right = root / "left", root / "right"
            left.mkdir()
            right.mkdir()
            names = ["B.fasta", "a.fasta", "c.fasta"]
            for directory, order in ((left, names), (right, reversed(names))):
                for name in order:
                    write_fasta(directory / name, [(name[0], "ACGT")])
            left_names = [path.name for path in batch.stable_paths(left, "*.fasta")]
            right_names = [path.name for path in batch.stable_paths(right, "*.fasta")]
            self.assertEqual(left_names, right_names)

            left_rows = batch.summarize_files(batch.stable_paths(left, "*.fasta"), "fasta", left / "summary.csv")
            right_rows = batch.summarize_files(batch.stable_paths(right, "*.fasta"), "fasta", right / "summary.csv")
            self.assertEqual([row["file"] for row in left_rows], [row["file"] for row in right_rows])
            self.assertEqual((left / "summary.csv").read_bytes(), (right / "summary.csv").read_bytes())

    def test_parallel_recipe_runs_under_spawn(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data = root / "data"
            data.mkdir()
            write_fasta(data / "b.fasta", [("b1", "AC"), ("b2", "ACGT")])
            write_fasta(data / "a.fasta", [("a1", "ACG")])
            launcher = root / "spawn_check.py"
            launcher.write_text(
                "from pathlib import Path\n"
                "import json, sys\n"
                f"sys.path.insert(0, {str(SKILL_ROOT / 'scripts')!r})\n"
                "from batch_process import process_files_parallel, stable_paths\n"
                "def main():\n"
                "    rows = process_files_parallel(stable_paths(Path(__file__).parent / 'data', '*.fasta'), workers=2, start_method='spawn')\n"
                "    print(json.dumps(rows, sort_keys=True))\n"
                "if __name__ == '__main__':\n"
                "    main()\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [sys.executable, "-B", str(launcher)],
                text=True,
                capture_output=True,
                timeout=20,
                check=True,
            )
            rows = json.loads(completed.stdout)
            self.assertEqual([row["file"] for row in rows], ["a.fasta", "b.fasta"])
            self.assertEqual(sum(row["count"] for row in rows), 3)
            self.assertEqual(sum(row["total_bp"] for row in rows), 9)
            fork_rows = batch.process_files_parallel(
                batch.stable_paths(data, "*.fasta"), workers=2, start_method="fork"
            )
            self.assertEqual(fork_rows, rows)

    def test_pyfastx_example_builds_and_reuses_gzip_index(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plain = root / "input.fasta"
            compressed = root / "input.fasta.gz"
            write_fasta(plain, [("alpha", "ACGT"), ("beta", "ACGTAC")])
            with plain.open("rb") as source, gzip.open(compressed, "wb") as target:
                target.write(source.read())

            first = pyfastx_example.inspect_fasta(compressed, "beta")
            second = pyfastx_example.inspect_fasta(compressed, "alpha")
            self.assertEqual(first["records"], 2)
            self.assertEqual(first["selected_length"], 6)
            self.assertEqual(second["selected_length"], 4)
            self.assertTrue(Path(first["index"]).is_file())
            completed = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(SKILL_ROOT / "scripts" / "pyfastx_index.py"),
                    str(compressed),
                    "--record-id",
                    "beta",
                ],
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertEqual(json.loads(completed.stdout)["selected_length"], 6)

    def test_standalone_demo_runs(self):
        completed = subprocess.run(
            [sys.executable, "-B", str(SKILL_ROOT / "scripts" / "batch_process.py")],
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertIn("Split sample0.fasta into 3 chunks", completed.stdout)
        self.assertIn("Indexed 15 records across 3 files", completed.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
