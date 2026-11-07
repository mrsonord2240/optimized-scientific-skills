#!/usr/bin/env python3
"""Run one analysis command with an owned process scope and staged publication.

The command is successful only when its root process exits zero, no process in
the owned scope remains, and a declared CSV passes a minimal structural check.
Failure receipts retain the staged files but never publish them. This is a
single-output transaction: the runner deliberately has no multi-output mode,
so it cannot expose partial multi-file publication or rollback semantics.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import shutil
import signal
import stat
import subprocess
import sys
import time
import uuid
from pathlib import Path


EXIT_VALIDATION = 65
EXIT_WORKER_FAILURE = 70
EXIT_LINGERING_CHILD = 75
EXIT_TIMEOUT = 124


class ScopeError(RuntimeError):
    """A process scope could not be created safely."""


class PosixScope:
    def __init__(self, argv: list[str]):
        self.process = subprocess.Popen(argv, start_new_session=True)
        self.root_pid = self.process.pid
        self.scope = f"process-group:{self.root_pid}"

    def wait(self, timeout: float) -> int:
        return self.process.wait(timeout=timeout)

    def active_after_root_exit(self) -> int:
        try:
            os.killpg(self.root_pid, 0)
        except ProcessLookupError:
            return 0
        except PermissionError as exc:
            raise ScopeError("cannot inspect owned process group") from exc
        return 1

    def cleanup(self) -> str:
        try:
            os.killpg(self.root_pid, signal.SIGTERM)
        except ProcessLookupError:
            return "none-needed"
        time.sleep(0.10)
        try:
            os.killpg(self.root_pid, signal.SIGKILL)
            return "SIGTERM-then-SIGKILL"
        except ProcessLookupError:
            return "SIGTERM"

    def close(self) -> None:
        return None


if os.name == "nt":
    import ctypes
    from ctypes import wintypes

    CREATE_SUSPENDED = 0x00000004
    CREATE_NEW_PROCESS_GROUP = 0x00000200
    CREATE_UNICODE_ENVIRONMENT = 0x00000400
    JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE = 0x00002000
    JobObjectBasicAccountingInformation = 1
    JobObjectExtendedLimitInformation = 9
    WAIT_OBJECT_0 = 0
    WAIT_TIMEOUT = 258
    INFINITE = 0xFFFFFFFF

    class STARTUPINFOW(ctypes.Structure):
        _fields_ = [
            ("cb", wintypes.DWORD), ("lpReserved", wintypes.LPWSTR),
            ("lpDesktop", wintypes.LPWSTR), ("lpTitle", wintypes.LPWSTR),
            ("dwX", wintypes.DWORD), ("dwY", wintypes.DWORD),
            ("dwXSize", wintypes.DWORD), ("dwYSize", wintypes.DWORD),
            ("dwXCountChars", wintypes.DWORD), ("dwYCountChars", wintypes.DWORD),
            ("dwFillAttribute", wintypes.DWORD), ("dwFlags", wintypes.DWORD),
            ("wShowWindow", wintypes.WORD), ("cbReserved2", wintypes.WORD),
            ("lpReserved2", ctypes.POINTER(ctypes.c_byte)),
            ("hStdInput", wintypes.HANDLE), ("hStdOutput", wintypes.HANDLE),
            ("hStdError", wintypes.HANDLE),
        ]

    class PROCESS_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("hProcess", wintypes.HANDLE), ("hThread", wintypes.HANDLE),
            ("dwProcessId", wintypes.DWORD), ("dwThreadId", wintypes.DWORD),
        ]

    class JOBOBJECT_BASIC_LIMIT_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("PerProcessUserTimeLimit", ctypes.c_longlong),
            ("PerJobUserTimeLimit", ctypes.c_longlong),
            ("LimitFlags", wintypes.DWORD),
            ("MinimumWorkingSetSize", ctypes.c_size_t),
            ("MaximumWorkingSetSize", ctypes.c_size_t),
            ("ActiveProcessLimit", wintypes.DWORD),
            ("Affinity", ctypes.c_size_t),
            ("PriorityClass", wintypes.DWORD),
            ("SchedulingClass", wintypes.DWORD),
        ]

    class IO_COUNTERS(ctypes.Structure):
        _fields_ = [(name, ctypes.c_ulonglong) for name in (
            "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
            "ReadTransferCount", "WriteTransferCount", "OtherTransferCount",
        )]

    class JOBOBJECT_EXTENDED_LIMIT_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("BasicLimitInformation", JOBOBJECT_BASIC_LIMIT_INFORMATION),
            ("IoInfo", IO_COUNTERS),
            ("ProcessMemoryLimit", ctypes.c_size_t),
            ("JobMemoryLimit", ctypes.c_size_t),
            ("PeakProcessMemoryUsed", ctypes.c_size_t),
            ("PeakJobMemoryUsed", ctypes.c_size_t),
        ]

    class JOBOBJECT_BASIC_ACCOUNTING_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("TotalUserTime", ctypes.c_longlong),
            ("TotalKernelTime", ctypes.c_longlong),
            ("ThisPeriodTotalUserTime", ctypes.c_longlong),
            ("ThisPeriodTotalKernelTime", ctypes.c_longlong),
            ("TotalPageFaultCount", wintypes.DWORD),
            ("TotalProcesses", wintypes.DWORD),
            ("ActiveProcesses", wintypes.DWORD),
            ("TotalTerminatedProcesses", wintypes.DWORD),
        ]

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
    kernel32.CreateJobObjectW.restype = wintypes.HANDLE
    kernel32.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
    kernel32.SetInformationJobObject.restype = wintypes.BOOL
    kernel32.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
    kernel32.AssignProcessToJobObject.restype = wintypes.BOOL
    kernel32.QueryInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD, ctypes.c_void_p]
    kernel32.QueryInformationJobObject.restype = wintypes.BOOL
    kernel32.CreateProcessW.argtypes = [
        wintypes.LPCWSTR, wintypes.LPWSTR, ctypes.c_void_p, ctypes.c_void_p,
        wintypes.BOOL, wintypes.DWORD, ctypes.c_void_p, wintypes.LPCWSTR,
        ctypes.POINTER(STARTUPINFOW), ctypes.POINTER(PROCESS_INFORMATION),
    ]
    kernel32.CreateProcessW.restype = wintypes.BOOL
    kernel32.ResumeThread.argtypes = [wintypes.HANDLE]
    kernel32.ResumeThread.restype = wintypes.DWORD
    kernel32.WaitForSingleObject.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    kernel32.WaitForSingleObject.restype = wintypes.DWORD
    kernel32.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
    kernel32.GetExitCodeProcess.restype = wintypes.BOOL
    kernel32.TerminateJobObject.argtypes = [wintypes.HANDLE, wintypes.UINT]
    kernel32.TerminateJobObject.restype = wintypes.BOOL
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel32.CloseHandle.restype = wintypes.BOOL

    def _wincheck(value: object, operation: str) -> None:
        if not value:
            raise ScopeError(f"{operation} failed (winerror {ctypes.get_last_error()})")

    class WindowsJobScope:
        def __init__(self, argv: list[str]):
            self.job = kernel32.CreateJobObjectW(None, None)
            _wincheck(self.job, "CreateJobObjectW")
            self.process = None
            self.thread = None
            try:
                limits = JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
                limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
                _wincheck(
                    kernel32.SetInformationJobObject(
                        self.job, JobObjectExtendedLimitInformation,
                        ctypes.byref(limits), ctypes.sizeof(limits),
                    ),
                    "SetInformationJobObject",
                )
                startup = STARTUPINFOW()
                startup.cb = ctypes.sizeof(startup)
                pi = PROCESS_INFORMATION()
                command_line = ctypes.create_unicode_buffer(subprocess.list2cmdline(argv))
                flags = CREATE_SUSPENDED | CREATE_NEW_PROCESS_GROUP | CREATE_UNICODE_ENVIRONMENT
                _wincheck(
                    kernel32.CreateProcessW(
                        None, command_line, None, None, False, flags, None, None,
                        ctypes.byref(startup), ctypes.byref(pi),
                    ),
                    "CreateProcessW",
                )
                self.process, self.thread, self.root_pid = pi.hProcess, pi.hThread, int(pi.dwProcessId)
                _wincheck(kernel32.AssignProcessToJobObject(self.job, self.process), "AssignProcessToJobObject")
                if kernel32.ResumeThread(self.thread) == 0xFFFFFFFF:
                    raise ScopeError(f"ResumeThread failed (winerror {ctypes.get_last_error()})")
                kernel32.CloseHandle(self.thread)
                self.thread = None
                self.scope = f"job-object:{self.root_pid}"
            except Exception:
                if self.process:
                    kernel32.TerminateJobObject(self.job, 1)
                self.close()
                raise

        def wait(self, timeout: float) -> int:
            result = kernel32.WaitForSingleObject(self.process, max(1, int(timeout * 1000)))
            if result == WAIT_TIMEOUT:
                raise subprocess.TimeoutExpired("owned-process", timeout)
            if result != WAIT_OBJECT_0:
                raise ScopeError(f"WaitForSingleObject failed (winerror {ctypes.get_last_error()})")
            code = wintypes.DWORD()
            _wincheck(kernel32.GetExitCodeProcess(self.process, ctypes.byref(code)), "GetExitCodeProcess")
            return int(code.value)

        def active_after_root_exit(self) -> int:
            info = JOBOBJECT_BASIC_ACCOUNTING_INFORMATION()
            _wincheck(
                kernel32.QueryInformationJobObject(
                    self.job, JobObjectBasicAccountingInformation,
                    ctypes.byref(info), ctypes.sizeof(info), None,
                ),
                "QueryInformationJobObject",
            )
            return int(info.ActiveProcesses)

        def cleanup(self) -> str:
            _wincheck(kernel32.TerminateJobObject(self.job, 1), "TerminateJobObject")
            return "TerminateJobObject"

        def close(self) -> None:
            if self.thread:
                kernel32.CloseHandle(self.thread)
                self.thread = None
            if self.process:
                kernel32.CloseHandle(self.process)
                self.process = None
            if self.job:
                kernel32.CloseHandle(self.job)
                self.job = None


def launch_scope(argv: list[str]):
    if os.name == "nt":
        return WindowsJobScope(argv)
    return PosixScope(argv)


def resolve_stage_path(stage: Path, relative: str) -> Path:
    relative_path = Path(relative)
    if relative_path.is_absolute():
        raise ValueError("staged paths must be relative to --stage-dir")
    # Resolve only for the containment check.  Returning the original
    # pathname lets the stable opener reject a final-component symlink rather
    # than silently treating its target as the declared artifact.
    candidate = stage / relative_path
    resolved = candidate.resolve()
    if resolved != stage and stage not in resolved.parents:
        raise ValueError("staged paths must remain beneath --stage-dir")
    return candidate


def _file_identity(info: os.stat_result) -> tuple[int, int]:
    return (info.st_dev, info.st_ino)


def read_stable_regular_bytes(path: Path) -> bytes:
    """Read one staged regular file without following a symlink or path race.

    A worker has exited before this is called, but the stage remains a path in
    a shared filesystem.  Hold a descriptor to the regular file through the
    read so validation and publication use the same bytes.
    """
    try:
        before = os.lstat(path)
    except OSError as exc:
        raise ValueError(f"missing CSV: {path}") from exc
    if not stat.S_ISREG(before.st_mode):
        raise ValueError(f"CSV must be a non-symlink regular file: {path}")

    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise ValueError(f"could not safely open CSV: {path}") from exc
    try:
        opened = os.fstat(descriptor)
        after = os.lstat(path)
        if (
            not stat.S_ISREG(opened.st_mode)
            or not stat.S_ISREG(after.st_mode)
            or _file_identity(before) != _file_identity(opened)
            or _file_identity(before) != _file_identity(after)
        ):
            raise ValueError(f"CSV changed or became non-regular while opening: {path}")
        with os.fdopen(descriptor, "rb") as handle:
            descriptor = None
            return handle.read()
    finally:
        if descriptor is not None:
            os.close(descriptor)


def validate_csv_bytes(payload: bytes, path: Path, required_columns: list[str], min_rows: int) -> dict[str, object]:
    if not payload:
        raise ValueError(f"missing or empty CSV: {path}")
    try:
        reader = csv.DictReader(io.StringIO(payload.decode("utf-8")))
        fields = reader.fieldnames or []
        missing = [column for column in required_columns if column not in fields]
        if missing:
            raise ValueError(f"CSV missing required columns: {', '.join(missing)}")
        row_count = sum(1 for _ in reader)
    except UnicodeDecodeError as exc:
        raise ValueError(f"CSV is not UTF-8: {path}") from exc
    if row_count < min_rows:
        raise ValueError(f"CSV has {row_count} rows; expected at least {min_rows}")
    return {"path": str(path), "columns": fields, "rows": row_count, "valid": True}


def write_receipt(path: Path, receipt: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def publish_verified_bytes_no_clobber(payload: bytes, destination: Path) -> None:
    """Atomically publish one new regular file from already validated bytes.

    Verified bytes are first fsynced to an invocation-owned sibling.  A
    create-only hard link then makes the complete file visible at the public
    name in one operation, refusing regular files and dangling symlinks alike.
    The temporary copy (not the stage file) is linked, so later stage mutation
    cannot change the public result.
    """
    if os.path.lexists(destination):
        raise FileExistsError(f"publish destination already exists: {destination}")
    descriptor: int | None = None
    temporary = destination.with_name(f".{destination.name}.{uuid.uuid4().hex}.tmp")
    temporary_identity: tuple[int, int] | None = None
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
        descriptor = os.open(temporary, flags, 0o600)
        temporary_identity = _file_identity(os.fstat(descriptor))
        with os.fdopen(descriptor, "wb") as handle:
            descriptor = None
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        current = os.lstat(temporary)
        if not stat.S_ISREG(current.st_mode) or _file_identity(current) != temporary_identity:
            raise OSError("temporary publication file changed or became non-regular")
        if os.path.lexists(destination):
            raise FileExistsError(f"publish destination already exists: {destination}")
        os.link(temporary, destination, follow_symlinks=False)
    except OSError as exc:
        if descriptor is not None:
            os.close(descriptor)
            descriptor = None
        raise OSError(f"refused to publish without an atomic no-clobber link: {exc}") from exc
    finally:
        # Remove only the regular temporary file created by this invocation;
        # never unlink a name replaced by another process.
        if temporary_identity is not None:
            try:
                current = os.lstat(temporary)
                if stat.S_ISREG(current.st_mode) and _file_identity(current) == temporary_identity:
                    temporary.unlink()
            except FileNotFoundError:
                pass


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", required=True, type=float, help="wall-clock seconds before owned cleanup")
    parser.add_argument("--grace", type=float, default=0.15, help="seconds to wait for owned descendants after root exit")
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--stage-dir", required=True, type=Path)
    parser.add_argument("--publish-source", required=True, help="CSV path relative to --stage-dir")
    parser.add_argument("--publish-dest", required=True, type=Path)
    parser.add_argument("--csv", required=True, help="CSV path relative to --stage-dir to validate")
    parser.add_argument("--require-columns", required=True, help="comma-separated CSV columns")
    parser.add_argument("--min-rows", type=int, default=1)
    parser.add_argument("command", nargs=argparse.REMAINDER, help="command to run, preceded by --")
    args = parser.parse_args(argv)
    if args.timeout <= 0 or args.grace < 0 or args.min_rows < 1:
        parser.error("timeout must be positive; grace non-negative; min-rows at least one")
    if args.command[:1] == ["--"]:
        args.command = args.command[1:]
    if not args.command:
        parser.error("provide a command after --")
    args.required_columns = [column.strip() for column in args.require_columns.split(",") if column.strip()]
    if not args.required_columns:
        parser.error("require at least one CSV column")
    return args


def run(args: argparse.Namespace) -> int:
    stage = args.stage_dir.resolve()
    receipt: dict[str, object] = {
        "command": args.command,
        "stage_dir": str(stage),
        "publish_source": args.publish_source,
        "publish_dest": str(args.publish_dest),
        "required_columns": args.required_columns,
        "min_rows": args.min_rows,
        "timeout_seconds": args.timeout,
        "grace_seconds": args.grace,
        "started_at": time.time(),
        "published": False,
    }
    scope = None
    try:
        # The declared stage must be created by this invocation.  Reusing an
        # existing directory would allow a successful worker that wrote
        # nothing to validate and publish a stale CSV left by an earlier run.
        # ``lexists`` also rejects a dangling symlink, for which ``exists`` is
        # false.  ``mkdir(exist_ok=False)`` closes the race after the
        # preflight.
        if os.path.lexists(args.stage_dir):
            raise ValueError(f"stage directory must not already exist: {args.stage_dir}")
        stage.parent.mkdir(parents=True, exist_ok=True)
        stage.mkdir(exist_ok=False)
        scope = launch_scope(args.command)
        receipt.update({"root_pid": scope.root_pid, "scope": scope.scope})
        try:
            worker_exit = scope.wait(args.timeout)
        except subprocess.TimeoutExpired:
            receipt.update({"outcome": "timeout", "cleanup": scope.cleanup()})
            return EXIT_TIMEOUT
        receipt["worker_exit"] = worker_exit
        if worker_exit != 0:
            receipt.update({"outcome": "worker-nonzero", "cleanup": scope.cleanup()})
            return EXIT_WORKER_FAILURE
        time.sleep(args.grace)
        active = scope.active_after_root_exit()
        receipt["active_processes_after_root_exit"] = active
        if active:
            receipt.update({"outcome": "lingering-owned-process", "cleanup": scope.cleanup()})
            return EXIT_LINGERING_CHILD
        csv_path = resolve_stage_path(stage, args.csv)
        source = resolve_stage_path(stage, args.publish_source)
        if csv_path != source:
            receipt.update({
                "outcome": "validation-failed",
                "validation_error": (
                    "--csv and --publish-source must resolve to the same staged file: "
                    f"{csv_path} != {source}"
                ),
            })
            return EXIT_VALIDATION
        try:
            payload = read_stable_regular_bytes(source)
            validated = validate_csv_bytes(payload, source, args.required_columns, args.min_rows)
        except (OSError, ValueError, csv.Error) as exc:
            receipt.update({"outcome": "validation-failed", "validation_error": str(exc)})
            return EXIT_VALIDATION
        args.publish_dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            publish_verified_bytes_no_clobber(payload, args.publish_dest)
        except OSError as exc:
            receipt.update({"outcome": "publish-refused", "publish_error": str(exc)})
            return EXIT_WORKER_FAILURE
        receipt.update({"outcome": "success", "validation": validated, "published": True})
        return 0
    except (OSError, ScopeError, ValueError) as exc:
        receipt.update({"outcome": "runner-error", "error": str(exc)})
        if scope is not None:
            try:
                receipt["cleanup"] = scope.cleanup()
            except ScopeError as cleanup_exc:
                receipt["cleanup_error"] = str(cleanup_exc)
        return EXIT_WORKER_FAILURE
    finally:
        receipt["finished_at"] = time.time()
        try:
            write_receipt(args.receipt, receipt)
        finally:
            if scope is not None:
                scope.close()


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(sys.argv[1:] if argv is None else argv))


if __name__ == "__main__":
    raise SystemExit(main())
