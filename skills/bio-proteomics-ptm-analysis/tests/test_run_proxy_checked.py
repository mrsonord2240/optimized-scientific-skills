"""Focused contract tests for the checked multi-artifact proxy runner."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
RUNNER = SCRIPTS / "run_proxy_checked.py"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("run_proxy_checked", RUNNER)
assert spec and spec.loader
proxy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(proxy)


SITES = (
    "Protein,Label,log2FC,adj.pvalue,pvalue_lfc,adj.pvalue_lfc,raw_ptm_log2FC,raw_ptm_adj.pvalue,adjustment_source,interpretation\n"
    "P1_S1,Treatment vs Control,1.5,0.01,0.02,0.03,1.6,0.02,co-enriched unmodified peptides; no paired global proteome,proxy-adjusted candidate; not a regulation call\n"
)
MODEL = "Protein,Label,log2FC,adj.pvalue\nP1_S1,Treatment vs Control,1.6,0.02\n"
QC = "Protein,Condition,BioReplicate,Raw.file,unmodified_rows,unique_unmodified_peptides\nP1,Control,1,C1,2,2\n"


def child_code(*, bad_hash: bool = False, exit_code: int = 0) -> str:
    return f'''from pathlib import Path
import csv, hashlib, os, sys
root = Path(sys.argv[1])
files = {{
    "proxy_adjusted_sites.csv": {SITES!r},
    "proxy_ptm_model.csv": {MODEL!r},
    "proxy_adjusted_qc.csv": {QC!r},
}}
for name, text in files.items():
    (root / name).write_text(text, encoding="utf-8")
rows = []
for name in ("proxy_adjusted_sites.csv", "proxy_ptm_model.csv", "proxy_adjusted_qc.csv"):
    path = root / name
    digest = hashlib.md5(path.read_bytes()).hexdigest()
    if {bad_hash!r} and name == "proxy_ptm_model.csv":
        digest = "0" * 32
    rows.append({{"schema_version": "1", "mode": "no-global-proxy", "complete": "TRUE", "artifact": name,
                 "bytes": str(path.stat().st_size), "md5": digest}})
with (root / "proxy_adjusted_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames={proxy.MANIFEST_COLUMNS!r})
    writer.writeheader(); writer.writerows(rows)
os._exit({exit_code})
'''


def write_valid_stage(root: Path) -> None:
    for name, text in ((proxy.ARTIFACTS[0], SITES), (proxy.ARTIFACTS[1], MODEL), (proxy.ARTIFACTS[2], QC)):
        (root / name).write_text(text, encoding="utf-8")
    rows = []
    for name in proxy.ARTIFACTS:
        path = root / name
        rows.append({"schema_version": "1", "mode": "no-global-proxy", "complete": "TRUE", "artifact": name,
                     "bytes": str(path.stat().st_size), "md5": hashlib.md5(path.read_bytes()).hexdigest()})
    with (root / proxy.MANIFEST).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=proxy.MANIFEST_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


class RunProxyCheckedTests(unittest.TestCase):
    def invoke(self, code: str, *, preexisting_publish: bool = False):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        base = Path(temporary.name)
        stage = base / "stage"
        publish_dir = base / "published"
        receipt = base / "receipt.json"
        if preexisting_publish:
            publish_dir.mkdir()
            (publish_dir / "old.csv").write_text("do not overwrite\n", encoding="utf-8")
        completed = subprocess.run(
            [
                sys.executable, str(RUNNER), "--timeout", "5", "--grace", "0.10",
                "--receipt", str(receipt), "--stage-dir", str(stage), "--publish-dir", str(publish_dir),
                "--", sys.executable, "-c", code, str(stage),
            ],
            capture_output=True, text=True, timeout=20,
        )
        return completed, json.loads(receipt.read_text(encoding="utf-8")), stage, publish_dir

    def test_success_publishes_verified_bundle(self):
        completed, receipt, _stage, publish_dir = self.invoke(child_code())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(receipt["outcome"], "success")
        self.assertTrue(receipt["published"])
        self.assertEqual({path.name for path in publish_dir.iterdir()}, {*proxy.ARTIFACTS, proxy.MANIFEST, proxy.OUTER_MARKER})
        marker = json.loads((publish_dir / proxy.OUTER_MARKER).read_text(encoding="utf-8"))
        self.assertEqual(marker["schema_version"], 1)
        self.assertEqual(marker["mode"], "no-global-proxy")
        self.assertTrue(marker["complete"])
        self.assertEqual(marker["worker_exit"], 0)
        self.assertEqual(marker["active_processes_after_root_exit"], 0)
        self.assertEqual(marker["internal_manifest"]["filename"], proxy.MANIFEST)
        self.assertEqual(receipt["validation"]["outer_marker"], proxy.OUTER_MARKER)

    def test_nonzero_child_leaves_no_public_outputs(self):
        completed, receipt, stage, publish_dir = self.invoke(child_code(exit_code=23))
        self.assertEqual(completed.returncode, 70)
        self.assertEqual(receipt["outcome"], "worker-nonzero")
        self.assertFalse(publish_dir.exists())
        self.assertFalse((publish_dir / proxy.OUTER_MARKER).exists())
        self.assertTrue((stage / proxy.MANIFEST).is_file())

    def test_timeout_leaves_no_public_outputs(self):
        completed, receipt, stage, publish_dir = self.invoke("import time; time.sleep(60)")
        self.assertEqual(completed.returncode, 124)
        self.assertEqual(receipt["outcome"], "timeout")
        self.assertFalse(publish_dir.exists())
        self.assertFalse((publish_dir / proxy.OUTER_MARKER).exists())
        self.assertTrue(stage.is_dir())

    def test_lingering_owned_child_leaves_no_public_outputs(self):
        code = (
            "import os, subprocess, sys, time; "
            "subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)']); "
            "time.sleep(.2); os._exit(0)"
        )
        completed, receipt, stage, publish_dir = self.invoke(code)
        self.assertEqual(completed.returncode, 75)
        self.assertEqual(receipt["outcome"], "lingering-owned-process")
        self.assertFalse(publish_dir.exists())
        self.assertFalse((publish_dir / proxy.OUTER_MARKER).exists())
        self.assertTrue(stage.is_dir())

    def test_invalid_manifest_hash_is_never_published(self):
        completed, receipt, stage, publish_dir = self.invoke(child_code(bad_hash=True))
        self.assertEqual(completed.returncode, 65)
        self.assertEqual(receipt["outcome"], "validation-failed")
        self.assertIn("hash or byte size", receipt["validation_error"])
        self.assertFalse(publish_dir.exists())
        self.assertFalse((publish_dir / proxy.OUTER_MARKER).exists())
        self.assertTrue((stage / proxy.MANIFEST).is_file())

    def test_existing_publish_directory_is_refused_without_launching_child(self):
        code = "from pathlib import Path; import sys; Path(sys.argv[1], 'worker-ran').write_text('bad', encoding='utf-8')"
        completed, receipt, stage, publish_dir = self.invoke(code, preexisting_publish=True)
        self.assertEqual(completed.returncode, 70)
        self.assertEqual(receipt["outcome"], "runner-error")
        self.assertIn("publish directory must not already exist", receipt["error"])
        self.assertEqual((publish_dir / "old.csv").read_text(encoding="utf-8"), "do not overwrite\n")
        self.assertFalse((publish_dir / proxy.OUTER_MARKER).exists())
        self.assertFalse(stage.exists())

    def test_outer_marker_is_written_last(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / "stage"
            stage.mkdir()
            for name, text in ((proxy.ARTIFACTS[0], SITES), (proxy.ARTIFACTS[1], MODEL), (proxy.ARTIFACTS[2], QC)):
                (stage / name).write_text(text, encoding="utf-8")
            rows = []
            for name in proxy.ARTIFACTS:
                path = stage / name
                rows.append({"schema_version": "1", "mode": "no-global-proxy", "complete": "TRUE", "artifact": name,
                             "bytes": str(path.stat().st_size), "md5": hashlib.md5(path.read_bytes()).hexdigest()})
            with (stage / proxy.MANIFEST).open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=proxy.MANIFEST_COLUMNS)
                writer.writeheader(); writer.writerows(rows)
            written: list[str] = []
            original_write = proxy.write_exclusive_bytes
            def recording_write(destination, payload):
                written.append(Path(destination).name)
                return original_write(destination, payload)
            with mock.patch.object(proxy, "write_exclusive_bytes", side_effect=recording_write):
                proxy.publish_proxy_bundle(stage, root / "published", worker_exit=0, active_processes=0)
            self.assertEqual(written, [proxy.MANIFEST, proxy.OUTER_MARKER])

    def test_outer_marker_cannot_be_created_without_clean_scope_state(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / "stage"
            publish_dir = root / "published"
            stage.mkdir()
            write_valid_stage(stage)
            with self.assertRaisesRegex(ValueError, "clean worker exit"):
                proxy.publish_proxy_bundle(stage, publish_dir, worker_exit=1, active_processes=0)
            self.assertFalse(publish_dir.exists())
            self.assertFalse((publish_dir / proxy.OUTER_MARKER).exists())

    def test_stage_mutation_after_success_does_not_change_public_bytes(self):
        completed, receipt, stage, publish_dir = self.invoke(child_code())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        public_model = publish_dir / proxy.ARTIFACTS[1]
        before = public_model.read_bytes()
        (stage / proxy.ARTIFACTS[1]).write_text("tampered after success\n", encoding="utf-8")
        self.assertEqual(public_model.read_bytes(), before)
        marker = json.loads((publish_dir / proxy.OUTER_MARKER).read_text(encoding="utf-8"))
        model = next(item for item in marker["artifacts"] if item["artifact"] == proxy.ARTIFACTS[1])
        self.assertEqual(hashlib.sha256(public_model.read_bytes()).hexdigest(), model["sha256"])
        self.assertTrue(receipt["published"])

    def test_stable_open_rejects_a_source_symlink_swap(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.csv"
            replacement = root / "replacement.csv"
            source.write_text("validated\n", encoding="utf-8")
            replacement.write_text("unvalidated\n", encoding="utf-8")
            real_open = os.open
            swapped = False

            def swap_before_open(path, flags, mode=0o777):
                nonlocal swapped
                if Path(path) == source and not swapped:
                    swapped = True
                    source.unlink()
                    source.symlink_to(replacement)
                return real_open(path, flags, mode)

            with mock.patch.object(proxy.os, "open", side_effect=swap_before_open):
                with self.assertRaisesRegex(ValueError, "changed or became non-regular|could not safely open"):
                    proxy.open_stable_regular(source)
            self.assertTrue(swapped)

    def test_copy_failure_surfaces_rollback_errors_and_never_writes_marker(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / "stage"
            stage.mkdir()
            for name, text in ((proxy.ARTIFACTS[0], SITES), (proxy.ARTIFACTS[1], MODEL), (proxy.ARTIFACTS[2], QC)):
                (stage / name).write_text(text, encoding="utf-8")
            rows = []
            for name in proxy.ARTIFACTS:
                path = stage / name
                rows.append({"schema_version": "1", "mode": "no-global-proxy", "complete": "TRUE", "artifact": name,
                             "bytes": str(path.stat().st_size), "md5": hashlib.md5(path.read_bytes()).hexdigest()})
            with (stage / proxy.MANIFEST).open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=proxy.MANIFEST_COLUMNS)
                writer.writeheader(); writer.writerows(rows)
            calls = 0
            original_copy = proxy.copy_verified_artifact
            def fail_second_copy(source, destination, expected_bytes, expected_md5):
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise OSError("synthetic copy failure")
                return original_copy(source, destination, expected_bytes, expected_md5)
            publish_dir = root / "published"
            locked_artifact = publish_dir / proxy.ARTIFACTS[0]
            real_unlink = Path.unlink
            def fail_first_rollback(path, *args, **kwargs):
                if path == locked_artifact:
                    raise OSError("synthetic locked output")
                return real_unlink(path, *args, **kwargs)
            with mock.patch.object(proxy, "copy_verified_artifact", side_effect=fail_second_copy), \
                 mock.patch.object(Path, "unlink", new=fail_first_rollback):
                with self.assertRaisesRegex(proxy.PublicationError, "rollback incomplete") as raised:
                    proxy.publish_proxy_bundle(stage, publish_dir, worker_exit=0, active_processes=0)
            self.assertTrue(raised.exception.rollback_errors)
            self.assertIn("synthetic locked output", " ".join(raised.exception.rollback_errors))
            self.assertTrue(locked_artifact.exists())
            self.assertFalse((publish_dir / proxy.OUTER_MARKER).exists())

    def test_runner_records_incomplete_rollback_in_its_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / "stage"
            publish_dir = root / "published"
            receipt_path = root / "receipt.json"

            class FinishedScope:
                root_pid = 12345
                scope = "test-scope"

                def wait(self, _timeout):
                    write_valid_stage(stage)
                    return 0

                def active_after_root_exit(self):
                    return 0

                def cleanup(self):
                    return "none-needed"

                def close(self):
                    return None

            args = proxy.parse_args([
                "--timeout", "5", "--receipt", str(receipt_path),
                "--stage-dir", str(stage), "--publish-dir", str(publish_dir), "--", "ignored-worker",
            ])
            calls = 0
            original_copy = proxy.copy_verified_artifact
            def fail_second_copy(source, destination, expected_bytes, expected_md5):
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise OSError("synthetic copy failure")
                return original_copy(source, destination, expected_bytes, expected_md5)
            locked_artifact = publish_dir / proxy.ARTIFACTS[0]
            real_unlink = Path.unlink
            def fail_first_rollback(path, *args, **kwargs):
                if path == locked_artifact:
                    raise OSError("synthetic locked output")
                return real_unlink(path, *args, **kwargs)
            with mock.patch.object(proxy, "launch_scope", return_value=FinishedScope()), \
                 mock.patch.object(proxy, "copy_verified_artifact", side_effect=fail_second_copy), \
                 mock.patch.object(Path, "unlink", new=fail_first_rollback):
                self.assertEqual(proxy.run(args), proxy.EXIT_WORKER_FAILURE)
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            self.assertEqual(receipt["outcome"], "publish-refused")
            self.assertIn("rollback incomplete", receipt["publish_error"])
            self.assertIn("synthetic locked output", " ".join(receipt["rollback_errors"]))
            self.assertFalse((publish_dir / proxy.OUTER_MARKER).exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
