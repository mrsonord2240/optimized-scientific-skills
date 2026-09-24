"""Focused red-capable tests for the stdlib-only checked-runner seam."""

from __future__ import annotations

import ctypes
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "run_checked.py"
spec = importlib.util.spec_from_file_location("run_checked", RUNNER)
assert spec and spec.loader
checked = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checked)


VALID_CSV = "site,score\nS1,1.0\n"
RUNNER_TIMEOUT_SECONDS = 5.0
INVOCATION_TIMEOUT_SECONDS = 20.0
CHILD_READY_TIMEOUT_SECONDS = 10.0
CHILD_CODE = (
    "from pathlib import Path; import os, sys, time; "
    "Path(sys.argv[1]).write_text(str(os.getpid()), encoding='utf-8'); time.sleep(60)"
)


def pid_is_alive(pid: int) -> bool:
    if os.name != "nt":
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        stat = Path(f"/proc/{pid}/stat")
        if stat.exists() and stat.read_text(encoding="utf-8").split()[2] == "Z":
            return False
        return True

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    process = kernel32.OpenProcess(0x00100000 | 0x1000, False, pid)  # SYNCHRONIZE | QUERY_LIMITED_INFORMATION
    if not process:
        return False
    try:
        code = ctypes.c_ulong()
        if not kernel32.GetExitCodeProcess(process, ctypes.byref(code)):
            return False
        return code.value == 259  # STILL_ACTIVE
    finally:
        kernel32.CloseHandle(process)


class RunCheckedTests(unittest.TestCase):
    def invoke(
        self,
        code: str,
        *,
        timeout: float = RUNNER_TIMEOUT_SECONDS,
        destination_kind: str | None = None,
        preseed_stage_csv: bool = False,
        csv_name: str = "result.csv",
        publish_source: str = "result.csv",
    ) -> tuple[subprocess.CompletedProcess[str], dict, Path, Path, Path]:
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        base = Path(temp.name)
        stage = base / "stage"
        destination = base / "published" / "result.csv"
        receipt = base / "receipt.json"
        if destination_kind == "file":
            destination.parent.mkdir(parents=True)
            destination.write_text("do not overwrite\n", encoding="utf-8")
        elif destination_kind == "dangling-symlink":
            destination.parent.mkdir(parents=True)
            try:
                destination.symlink_to(base / "missing-target.csv")
            except (NotImplementedError, OSError) as exc:
                self.skipTest(f"symlink fixture unavailable: {exc}")
        if preseed_stage_csv:
            stage.mkdir(parents=True)
            (stage / "result.csv").write_text(VALID_CSV, encoding="utf-8")
        completed = subprocess.run(
            [
                sys.executable, str(RUNNER), "--timeout", str(timeout), "--grace", "0.10",
                "--receipt", str(receipt), "--stage-dir", str(stage),
                "--publish-source", publish_source, "--publish-dest", str(destination),
                "--csv", csv_name, "--require-columns", "site,score", "--min-rows", "1",
                "--", sys.executable, "-c", code, str(stage),
            ],
            capture_output=True,
            text=True,
            timeout=INVOCATION_TIMEOUT_SECONDS,
        )
        return completed, json.loads(receipt.read_text(encoding="utf-8")), stage, destination, base

    def assert_pid_stops(self, pid: int) -> None:
        deadline = time.monotonic() + 2.0
        while time.monotonic() < deadline and pid_is_alive(pid):
            time.sleep(0.02)
        self.assertFalse(pid_is_alive(pid), f"owned child PID {pid} survived cleanup")

    def test_success_publishes_only_after_valid_clean_completion(self) -> None:
        code = "from pathlib import Path; import sys; p=Path(sys.argv[1]); p.mkdir(exist_ok=True); (p/'result.csv').write_text(" + repr(VALID_CSV) + ")"
        completed, receipt, _stage, destination, _base = self.invoke(code)
        self.assertEqual(completed.returncode, 0, f"{completed.stderr}\n{receipt}")
        self.assertEqual(receipt["outcome"], "success")
        self.assertTrue(receipt["published"])
        self.assertEqual(receipt["publish_source"], "result.csv")
        self.assertEqual(receipt["publish_dest"], str(destination))
        self.assertEqual(receipt["required_columns"], ["site", "score"])
        self.assertEqual(receipt["min_rows"], 1)
        self.assertEqual(receipt["timeout_seconds"], RUNNER_TIMEOUT_SECONDS)
        self.assertEqual(receipt["grace_seconds"], 0.10)
        self.assertEqual(destination.read_text(encoding="utf-8"), VALID_CSV)

    def test_publication_copies_verified_bytes_instead_of_aliasing_stage(self) -> None:
        code = "from pathlib import Path; import sys; p=Path(sys.argv[1]); p.mkdir(exist_ok=True); (p/'result.csv').write_text(" + repr(VALID_CSV) + ")"
        completed, receipt, stage, destination, _base = self.invoke(code)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(receipt["published"])
        (stage / "result.csv").write_text("site,score\nS2,2.0\n", encoding="utf-8")
        self.assertEqual(destination.read_text(encoding="utf-8"), VALID_CSV)
        self.assertEqual(list(destination.parent.glob(".result.csv.*.tmp")), [])

    def test_output_then_139_is_not_published(self) -> None:
        code = "from pathlib import Path; import os, sys; p=Path(sys.argv[1]); p.mkdir(exist_ok=True); (p/'result.csv').write_text(" + repr(VALID_CSV) + "); os._exit(139)"
        completed, receipt, stage, destination, _base = self.invoke(code)
        self.assertEqual(completed.returncode, 70)
        self.assertEqual(receipt["outcome"], "worker-nonzero")
        self.assertEqual(receipt["worker_exit"], 139)
        self.assertEqual(receipt["publish_dest"], str(destination))
        self.assertEqual(receipt["required_columns"], ["site", "score"])
        self.assertFalse(destination.exists())
        self.assertTrue((stage / "result.csv").exists())

    def test_preexisting_stage_cannot_publish_stale_csv_and_gets_a_receipt(self) -> None:
        # A stage is a transaction boundary owned by one invocation.  If a
        # valid-looking CSV is already there, reject before the worker starts;
        # otherwise a zero-exit worker that wrote nothing could publish stale
        # data from a previous run.
        code = (
            "from pathlib import Path; import sys; "
            "Path(sys.argv[1], 'worker-ran').write_text('unexpected', encoding='utf-8')"
        )
        completed, receipt, stage, destination, _base = self.invoke(
            code,
            preseed_stage_csv=True,
        )
        self.assertEqual(completed.returncode, 70)
        self.assertEqual(receipt["outcome"], "runner-error")
        self.assertIn("stage directory must not already exist", receipt["error"])
        self.assertFalse(destination.exists())
        self.assertEqual((stage / "result.csv").read_text(encoding="utf-8"), VALID_CSV)
        self.assertFalse((stage / "worker-ran").exists())

    def test_timeout_cleans_only_the_owned_child_tree(self) -> None:
        code = (
            "from pathlib import Path; import subprocess, sys, time; p=Path(sys.argv[1]); p.mkdir(exist_ok=True); "
            "subprocess.Popen([sys.executable, '-c', " + repr(CHILD_CODE) + ", str(p/'child.pid')]); "
            "deadline=time.monotonic()+" + repr(CHILD_READY_TIMEOUT_SECONDS) + "; "
            "exec(\"while not (p/'child.pid').exists() and time.monotonic() < deadline:\\n time.sleep(.01)\"); "
            "assert (p/'child.pid').exists(), 'child did not signal readiness'; "
            "(p/'parent.ready').write_text('ready', encoding='utf-8'); time.sleep(60)"
        )
        completed, receipt, stage, destination, _base = self.invoke(code)
        self.assertEqual(completed.returncode, 124)
        self.assertEqual(receipt["outcome"], "timeout")
        self.assertFalse(destination.exists())
        self.assertTrue((stage / "parent.ready").is_file())
        self.assert_pid_stops(int((stage / "child.pid").read_text(encoding="utf-8")))

    def test_lingering_child_after_zero_exit_is_not_published(self) -> None:
        code = (
            "from pathlib import Path; import os, subprocess, sys, time; p=Path(sys.argv[1]); p.mkdir(exist_ok=True); "
            "(p/'result.csv').write_text(" + repr(VALID_CSV) + "); "
            "subprocess.Popen([sys.executable, '-c', " + repr(CHILD_CODE) + ", str(p/'child.pid')]); "
            "deadline=time.monotonic()+" + repr(CHILD_READY_TIMEOUT_SECONDS) + "; "
            "exec(\"while not (p/'child.pid').exists() and time.monotonic() < deadline:\\n time.sleep(.01)\"); "
            "assert (p/'child.pid').exists(), 'child did not signal readiness'; "
            "(p/'parent.ready').write_text('ready', encoding='utf-8'); os._exit(0)"
        )
        completed, receipt, stage, destination, _base = self.invoke(code)
        self.assertEqual(completed.returncode, 75)
        self.assertEqual(receipt["outcome"], "lingering-owned-process")
        self.assertFalse(destination.exists())
        self.assertTrue((stage / "parent.ready").is_file())
        self.assert_pid_stops(int((stage / "child.pid").read_text(encoding="utf-8")))

    def test_malformed_csv_is_not_published(self) -> None:
        code = "from pathlib import Path; import sys; p=Path(sys.argv[1]); p.mkdir(exist_ok=True); (p/'result.csv').write_text('wrong\\n1\\n')"
        completed, receipt, stage, destination, _base = self.invoke(code)
        self.assertEqual(completed.returncode, 65)
        self.assertEqual(receipt["outcome"], "validation-failed")
        self.assertFalse(destination.exists())
        self.assertTrue((stage / "result.csv").exists())

    def test_publish_source_must_be_the_validated_staged_file(self) -> None:
        code = (
            "from pathlib import Path; import sys; p=Path(sys.argv[1]); p.mkdir(exist_ok=True); "
            "(p/'validated.csv').write_text(" + repr(VALID_CSV) + "); "
            "(p/'published.csv').write_text(" + repr(VALID_CSV) + ")"
        )
        completed, receipt, stage, destination, _base = self.invoke(
            code,
            csv_name="validated.csv",
            publish_source="published.csv",
        )
        self.assertEqual(completed.returncode, 65)
        self.assertEqual(receipt["outcome"], "validation-failed")
        self.assertIn("--csv and --publish-source must resolve to the same staged file", receipt["validation_error"])
        self.assertFalse(destination.exists())
        self.assertTrue((stage / "validated.csv").exists())
        self.assertTrue((stage / "published.csv").exists())

    def test_existing_destination_is_never_clobbered(self) -> None:
        code = "from pathlib import Path; import sys; p=Path(sys.argv[1]); p.mkdir(exist_ok=True); (p/'result.csv').write_text(" + repr(VALID_CSV) + ")"
        completed, receipt, stage, destination, _base = self.invoke(code, destination_kind="file")
        self.assertEqual(completed.returncode, 70)
        self.assertEqual(receipt["outcome"], "publish-refused")
        self.assertEqual(destination.read_text(encoding="utf-8"), "do not overwrite\n")
        self.assertTrue((stage / "result.csv").exists())

    def test_dangling_symlink_destination_is_never_replaced(self) -> None:
        code = "from pathlib import Path; import sys; p=Path(sys.argv[1]); p.mkdir(exist_ok=True); (p/'result.csv').write_text(" + repr(VALID_CSV) + ")"
        completed, receipt, stage, destination, _base = self.invoke(code, destination_kind="dangling-symlink")
        self.assertEqual(completed.returncode, 70)
        self.assertEqual(receipt["outcome"], "publish-refused")
        self.assertTrue(destination.is_symlink())
        self.assertFalse(destination.exists())
        self.assertTrue((stage / "result.csv").exists())

    def test_staged_source_symlink_is_refused(self) -> None:
        code = (
            "from pathlib import Path; import sys; p=Path(sys.argv[1]); p.mkdir(exist_ok=True); "
            "(p/'payload.csv').write_text(" + repr(VALID_CSV) + "); "
            "(p/'result.csv').symlink_to(p/'payload.csv')"
        )
        completed, receipt, stage, destination, _base = self.invoke(code)
        if not (stage / "result.csv").is_symlink():
            self.skipTest("symlink fixture unavailable")
        self.assertEqual(completed.returncode, 65)
        self.assertEqual(receipt["outcome"], "validation-failed")
        self.assertIn("non-symlink regular file", receipt["validation_error"])
        self.assertTrue((stage / "result.csv").is_symlink())
        self.assertFalse(destination.exists())

    def test_stable_open_rejects_stage_path_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "result.csv"
            replacement = root / "replacement.csv"
            source.write_text(VALID_CSV, encoding="utf-8")
            replacement.write_text("site,score\nS2,2.0\n", encoding="utf-8")
            real_open = os.open
            swapped = False

            def replace_before_open(path, flags, mode=0o777):
                nonlocal swapped
                if Path(path) == source and not swapped:
                    swapped = True
                    source.unlink()
                    os.replace(replacement, source)
                return real_open(path, flags, mode)

            with mock.patch.object(checked.os, "open", side_effect=replace_before_open):
                with self.assertRaisesRegex(ValueError, "changed or became non-regular|could not safely open"):
                    checked.read_stable_regular_bytes(source)
            self.assertTrue(swapped)


if __name__ == "__main__":
    unittest.main(verbosity=2)
