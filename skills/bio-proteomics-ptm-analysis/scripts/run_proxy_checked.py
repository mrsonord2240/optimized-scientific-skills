#!/usr/bin/env python3
"""Run a no-global PTM proxy child and publish a checked, marker-last bundle.

The worker's stage is never aliased into the public result.  After clean owned
process completion, this runner opens stable regular stage files, copies and
checks them into an invocation-owned directory, and writes its own completion
marker last.  A worker-written manifest alone is therefore diagnostic output,
not evidence of checked publication.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import stat
import subprocess
import sys
import time
import uuid
from pathlib import Path

from run_checked import (
    EXIT_LINGERING_CHILD,
    EXIT_TIMEOUT,
    EXIT_VALIDATION,
    EXIT_WORKER_FAILURE,
    ScopeError,
    launch_scope,
)


ARTIFACTS = (
    "proxy_adjusted_sites.csv",
    "proxy_ptm_model.csv",
    "proxy_adjusted_qc.csv",
)
MANIFEST = "proxy_adjusted_manifest.csv"
OUTER_MARKER = "proxy_checked_complete.json"
MANIFEST_COLUMNS = ("schema_version", "mode", "complete", "artifact", "bytes", "md5")
SITE_COLUMNS = (
    "Protein", "Label", "log2FC", "adj.pvalue", "pvalue_lfc", "adj.pvalue_lfc",
    "raw_ptm_log2FC", "raw_ptm_adj.pvalue", "adjustment_source", "interpretation",
)
MODEL_COLUMNS = ("Protein", "Label", "log2FC", "adj.pvalue")
QC_COLUMNS = ("Protein", "Condition", "BioReplicate", "Raw.file", "unmodified_rows", "unique_unmodified_peptides")
ADJUSTMENT_SOURCE = "co-enriched unmodified peptides; no paired global proteome"
INTERPRETATION = "proxy-adjusted candidate; not a regulation call"


class PublicationError(OSError):
    """A publish attempt failed, optionally with rollback failures to report."""

    def __init__(self, message: str, rollback_errors: list[str] | None = None):
        super().__init__(message)
        self.rollback_errors = rollback_errors or []


def _file_identity(info: os.stat_result) -> tuple[int, int]:
    return (info.st_dev, info.st_ino)


def open_stable_regular(path: Path) -> int:
    """Open a non-symlink regular file and detect pathname replacement races.

    ``O_NOFOLLOW`` closes the lstat/open race where the platform supports it.
    Windows does not expose that flag in Python, so identity checks before and
    after open fail closed if a pathname becomes a link or another file.
    """
    try:
        before = os.lstat(path)
    except OSError as exc:
        raise ValueError(f"missing proxy file: {path.name}") from exc
    if not stat.S_ISREG(before.st_mode):
        raise ValueError(f"proxy file must be a non-symlink regular file: {path.name}")

    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise ValueError(f"could not safely open proxy file: {path.name}") from exc
    try:
        opened = os.fstat(descriptor)
        after = os.lstat(path)
        if (
            not stat.S_ISREG(opened.st_mode)
            or not stat.S_ISREG(after.st_mode)
            or _file_identity(before) != _file_identity(opened)
            or _file_identity(before) != _file_identity(after)
        ):
            raise ValueError(f"proxy file changed or became non-regular while opening: {path.name}")
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def read_stable_bytes(path: Path) -> bytes:
    descriptor = open_stable_regular(path)
    try:
        with os.fdopen(descriptor, "rb") as handle:
            return handle.read()
    except OSError as exc:
        raise ValueError(f"could not read proxy file: {path.name}") from exc


def read_csv_rows(path: Path, required_columns: tuple[str, ...]) -> list[dict[str, str]]:
    raw = read_stable_bytes(path)
    if not raw:
        raise ValueError(f"missing, empty, or non-regular proxy artifact: {path.name}")
    try:
        reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
        fields = tuple(reader.fieldnames or ())
        missing = [column for column in required_columns if column not in fields]
        if missing:
            raise ValueError(f"{path.name} missing required columns: {', '.join(missing)}")
        rows = list(reader)
    except UnicodeDecodeError as exc:
        raise ValueError(f"{path.name} is not UTF-8 CSV") from exc
    if not rows:
        raise ValueError(f"{path.name} has no data rows")
    return rows


def read_proxy_manifest(stage: Path) -> tuple[list[dict[str, str]], bytes]:
    """Read and validate the small worker manifest from a stable file handle."""
    manifest_path = stage / MANIFEST
    manifest_bytes = read_stable_bytes(manifest_path)
    if not manifest_bytes:
        raise ValueError(f"missing, empty, or non-regular proxy manifest: {MANIFEST}")
    try:
        reader = csv.DictReader(io.StringIO(manifest_bytes.decode("utf-8")))
        if tuple(reader.fieldnames or ()) != MANIFEST_COLUMNS:
            raise ValueError("proxy manifest has an unexpected schema")
        manifest_rows = list(reader)
    except UnicodeDecodeError as exc:
        raise ValueError("proxy manifest is not UTF-8 CSV") from exc
    if len(manifest_rows) != len(ARTIFACTS):
        raise ValueError("proxy manifest must contain exactly three artifact rows")

    manifest_names = tuple(row.get("artifact", "") for row in manifest_rows)
    if manifest_names != ARTIFACTS:
        raise ValueError("proxy manifest artifact basenames must be exactly the expected three in order")

    for index, name in enumerate(ARTIFACTS):
        row = manifest_rows[index]
        if row.get("schema_version") != "1" or row.get("mode") != "no-global-proxy" or row.get("complete") != "TRUE":
            raise ValueError(f"proxy manifest row for {name} has an invalid schema/version/mode/complete value")
        try:
            declared_bytes = int(row["bytes"])
        except (KeyError, ValueError) as exc:
            raise ValueError(f"proxy manifest bytes for {name} is invalid") from exc
        if declared_bytes < 0:
            raise ValueError(f"proxy manifest bytes for {name} is invalid")
        declared_md5 = row.get("md5", "")
        if len(declared_md5) != 32 or any(character not in "0123456789abcdefABCDEF" for character in declared_md5):
            raise ValueError(f"proxy manifest MD5 for {name} is invalid")
    return manifest_rows, manifest_bytes


def _declared_artifact_metadata(manifest_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    """Normalize the already schema-checked worker manifest for publication."""
    artifact_metadata: list[dict[str, object]] = []
    # Indexing keeps this stdlib helper compatible with the repository's
    # Python 3.9+ support boundary (``zip(..., strict=True)`` needs 3.10).
    for index, name in enumerate(ARTIFACTS):
        row = manifest_rows[index]
        if row.get("schema_version") != "1" or row.get("mode") != "no-global-proxy" or row.get("complete") != "TRUE":
            raise ValueError(f"proxy manifest row for {name} has an invalid schema/version/mode/complete value")
        try:
            declared_bytes = int(row["bytes"])
        except (KeyError, ValueError) as exc:
            raise ValueError(f"proxy manifest bytes for {name} is invalid") from exc
        declared_md5 = row.get("md5", "")
        if len(declared_md5) != 32 or any(character not in "0123456789abcdefABCDEF" for character in declared_md5):
            raise ValueError(f"proxy manifest MD5 for {name} is invalid")
        artifact_metadata.append({"artifact": name, "bytes": declared_bytes, "md5": declared_md5.casefold()})
    return artifact_metadata


def validate_proxy_bundle(stage: Path) -> dict[str, object]:
    """Validate a stage in place for diagnostics; publishing rechecks during copy."""
    manifest_rows, _manifest_bytes = read_proxy_manifest(stage)
    artifact_metadata = _declared_artifact_metadata(manifest_rows)
    for metadata in artifact_metadata:
        source = stage / str(metadata["artifact"])
        descriptor = open_stable_regular(source)
        digest = hashlib.md5()
        size = 0
        with os.fdopen(descriptor, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
        if size != metadata["bytes"] or digest.hexdigest() != metadata["md5"]:
            raise ValueError(f"proxy manifest hash or byte size does not match {source.name}")

    site_rows = read_csv_rows(stage / ARTIFACTS[0], SITE_COLUMNS)
    for row in site_rows:
        if row["adjustment_source"] != ADJUSTMENT_SOURCE or row["interpretation"] != INTERPRETATION:
            raise ValueError("proxy-adjusted sites are missing the required non-regulatory labels")
    read_csv_rows(stage / ARTIFACTS[1], MODEL_COLUMNS)
    read_csv_rows(stage / ARTIFACTS[2], QC_COLUMNS)
    return {"manifest": MANIFEST, "artifacts": artifact_metadata, "site_rows": len(site_rows), "valid": True}


def copy_verified_artifact(source: Path, destination: Path, expected_bytes: int, expected_md5: str) -> dict[str, object]:
    """Copy a stable stage file to a create-only destination while hashing it."""
    source_descriptor = open_stable_regular(source)
    destination_descriptor: int | None = None
    digest_md5 = hashlib.md5()
    digest_sha256 = hashlib.sha256()
    copied = 0
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
        destination_descriptor = os.open(destination, flags, 0o600)
        with os.fdopen(source_descriptor, "rb") as source_handle:
            source_descriptor = None
            with os.fdopen(destination_descriptor, "wb") as destination_handle:
                destination_descriptor = None
                for chunk in iter(lambda: source_handle.read(1024 * 1024), b""):
                    copied += len(chunk)
                    digest_md5.update(chunk)
                    digest_sha256.update(chunk)
                    destination_handle.write(chunk)
                destination_handle.flush()
                os.fsync(destination_handle.fileno())
        if copied != expected_bytes or digest_md5.hexdigest() != expected_md5.casefold():
            raise ValueError(f"proxy manifest hash or byte size does not match {source.name}")
        return {"artifact": source.name, "bytes": copied, "md5": digest_md5.hexdigest(), "sha256": digest_sha256.hexdigest()}
    finally:
        if source_descriptor is not None:
            os.close(source_descriptor)
        if destination_descriptor is not None:
            os.close(destination_descriptor)


def write_exclusive_bytes(destination: Path, payload: bytes) -> dict[str, object]:
    """Write a new regular file without following or replacing a name."""
    descriptor: int | None = None
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
        descriptor = os.open(destination, flags, 0o600)
        with os.fdopen(descriptor, "wb") as handle:
            descriptor = None
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        return {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
    finally:
        if descriptor is not None:
            os.close(descriptor)


def _rollback_publish(publish_dir: Path, published: list[Path]) -> list[str]:
    errors: list[str] = []
    for destination in reversed(published):
        try:
            destination.unlink()
        except FileNotFoundError:
            pass
        except OSError as exc:
            errors.append(f"could not remove {destination.name}: {exc}")
    try:
        publish_dir.rmdir()
    except FileNotFoundError:
        pass
    except OSError as exc:
        errors.append(f"could not remove publish directory: {exc}")
    return errors


def publish_proxy_bundle(stage: Path, publish_dir: Path, *, worker_exit: int, active_processes: int) -> dict[str, object]:
    """Copy a checked bundle into a fresh directory and write the outer marker last."""
    if worker_exit != 0 or active_processes != 0:
        raise ValueError("refusing proxy publication before clean worker exit and zero owned processes")
    if os.path.lexists(publish_dir):
        raise FileExistsError(f"proxy publish directory already exists: {publish_dir}")
    publish_dir.parent.mkdir(parents=True, exist_ok=True)
    try:
        publish_dir.mkdir()
    except FileExistsError as exc:
        raise FileExistsError(f"proxy publish directory already exists: {publish_dir}") from exc

    published: list[Path] = []
    try:
        manifest_rows, manifest_bytes = read_proxy_manifest(stage)
        declared_artifacts = _declared_artifact_metadata(manifest_rows)
        copied_artifacts: list[dict[str, object]] = []
        for metadata in declared_artifacts:
            name = str(metadata["artifact"])
            source = stage / name
            destination = publish_dir / name
            published.append(destination)
            copied_artifacts.append(copy_verified_artifact(source, destination, int(metadata["bytes"]), str(metadata["md5"])))

        site_rows = read_csv_rows(publish_dir / ARTIFACTS[0], SITE_COLUMNS)
        for row in site_rows:
            if row["adjustment_source"] != ADJUSTMENT_SOURCE or row["interpretation"] != INTERPRETATION:
                raise ValueError("proxy-adjusted sites are missing the required non-regulatory labels")
        read_csv_rows(publish_dir / ARTIFACTS[1], MODEL_COLUMNS)
        read_csv_rows(publish_dir / ARTIFACTS[2], QC_COLUMNS)

        manifest_destination = publish_dir / MANIFEST
        published.append(manifest_destination)
        manifest_copy = write_exclusive_bytes(manifest_destination, manifest_bytes)
        marker_payload = {
            "schema_version": 1,
            "mode": "no-global-proxy",
            "complete": True,
            "worker_exit": worker_exit,
            "active_processes_after_root_exit": active_processes,
            "internal_manifest": {"filename": MANIFEST, **manifest_copy},
            "artifacts": copied_artifacts,
        }
        marker_destination = publish_dir / OUTER_MARKER
        published.append(marker_destination)
        marker_copy = write_exclusive_bytes(
            marker_destination,
            (json.dumps(marker_payload, indent=2, sort_keys=True) + "\n").encode("utf-8"),
        )
        return {
            "manifest": MANIFEST,
            "outer_marker": OUTER_MARKER,
            "artifacts": copied_artifacts,
            "site_rows": len(site_rows),
            "marker": marker_copy,
            "valid": True,
        }
    except Exception as exc:
        rollback_errors = _rollback_publish(publish_dir, published)
        message = f"proxy publication failed: {exc}"
        if rollback_errors:
            message += "; rollback incomplete: " + "; ".join(rollback_errors)
        raise PublicationError(message, rollback_errors) from exc


def write_receipt(path: Path, receipt: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", required=True, type=float)
    parser.add_argument("--grace", type=float, default=0.15)
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--stage-dir", required=True, type=Path)
    parser.add_argument("--publish-dir", required=True, type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER, help="command to run, preceded by --")
    args = parser.parse_args(argv)
    if args.timeout <= 0 or args.grace < 0:
        parser.error("timeout must be positive and grace non-negative")
    if args.command[:1] == ["--"]:
        args.command = args.command[1:]
    if not args.command:
        parser.error("provide a command after --")
    return args


def run(args: argparse.Namespace) -> int:
    stage = args.stage_dir.resolve()
    publish_dir = args.publish_dir.resolve()
    receipt: dict[str, object] = {
        "command": args.command,
        "stage_dir": str(stage),
        "publish_dir": str(publish_dir),
        "timeout_seconds": args.timeout,
        "grace_seconds": args.grace,
        "started_at": time.time(),
        "published": False,
    }
    scope = None
    try:
        if os.path.lexists(args.stage_dir):
            raise ValueError(f"stage directory must not already exist: {args.stage_dir}")
        if os.path.lexists(args.publish_dir):
            raise ValueError(f"proxy publish directory must not already exist: {args.publish_dir}")
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
        try:
            validation = validate_proxy_bundle(stage)
        except (OSError, ValueError, csv.Error) as exc:
            receipt.update({"outcome": "validation-failed", "validation_error": str(exc)})
            return EXIT_VALIDATION
        try:
            validation = publish_proxy_bundle(
                stage,
                publish_dir,
                worker_exit=worker_exit,
                active_processes=active,
            )
        except PublicationError as exc:
            receipt.update({"outcome": "publish-refused", "publish_error": str(exc)})
            if exc.rollback_errors:
                receipt["rollback_errors"] = exc.rollback_errors
            return EXIT_WORKER_FAILURE
        receipt.update({"outcome": "success", "validation": validation, "published": True})
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
