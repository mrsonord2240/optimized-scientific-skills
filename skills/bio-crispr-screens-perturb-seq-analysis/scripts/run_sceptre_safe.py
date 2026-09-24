#!/usr/bin/env python3
"""Run sceptre in a child process and accept only a validated result table.

Some Windows R/sceptre stacks fault during native teardown after writing a table.
The manifest makes that worker status visible while keeping the caller contract tied
to a parseable, complete result rather than to the child process's teardown status.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path


def validate_result(path: Path) -> tuple[bool, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return False, "result table was not written"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "p_value" not in reader.fieldnames:
            return False, "result table has no p_value column"
        rows = list(reader)
    if not rows:
        return False, "result table has no data rows"
    for row in rows:
        value = row.get("p_value", "")
        if value:
            try:
                p_value = float(value)
            except ValueError:
                return False, f"non-numeric p_value: {value!r}"
            if not 0.0 <= p_value <= 1.0:
                return False, f"p_value outside [0, 1]: {p_value}"
    return True, f"validated {len(rows)} rows"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", help="input .rds (omit with --example)")
    parser.add_argument("out", help="output TSV")
    parser.add_argument("--example", action="store_true", help="run sceptre bundled low-MOI example")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable or project R wrapper")
    args = parser.parse_args()
    if args.example == (args.input is not None):
        parser.error("use exactly one of --example or input.rds")

    script = Path(__file__).with_name("run_sceptre.R")
    out = Path(args.out)
    rscript_command = shlex.split(args.rscript)
    # Project R launchers are often POSIX shell files that set R_LIBS_USER and PATH.
    # Invoke those through bash instead of asking Windows CreateProcess to execute a .sh file.
    shell_wrapper = len(rscript_command) == 1 and rscript_command[0].lower().endswith(".sh")
    if shell_wrapper:
        # MSYS may convert /f/path passed to native Python into F:/path. Convert it back
        # because the nested bash invocation needs a POSIX path for the wrapper itself.
        match = re.match(r"^([A-Za-z]):/(.*)$", rscript_command[0])
        if match:
            rscript_command[0] = f"/{match.group(1).lower()}/{match.group(2)}"
        git_bash = Path(os.environ.get("GIT_BASH", r"C:\Program Files\Git\bin\bash.exe"))
        rscript_command.insert(0, str(git_bash) if git_bash.is_file() else "bash")
    command = [*rscript_command, str(script)]
    command.extend(["--example", str(out)] if args.example else [args.input, str(out)])
    try:
        environment = os.environ.copy()
        if shell_wrapper:
            # Preserve /f/... POSIX paths when native Python launches Git Bash.
            environment["MSYS_NO_PATHCONV"] = "1"
        completed = subprocess.run(command, text=True, env=environment)
    except OSError as error:
        completed = subprocess.CompletedProcess(command, returncode=127)
        launch_error = str(error)
    else:
        launch_error = None
    valid, detail = validate_result(out)
    manifest = {
        "command": command,
        "worker_returncode": completed.returncode,
        "result_validation": detail,
        "accepted": valid,
        "launch_error": launch_error,
        "warning": (
            "worker returned non-zero after a validated result; inspect this manifest and pin a stable R/sceptre stack"
            if valid and completed.returncode else None
        ),
    }
    manifest_path = out.with_suffix(out.suffix + ".run.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if not valid:
        print(f"SCEPTRE result rejected: {detail}; worker exit={completed.returncode}", file=sys.stderr)
        return completed.returncode or 1
    if completed.returncode:
        print(f"SCEPTRE worker exit={completed.returncode}, but {detail}; details: {manifest_path}", file=sys.stderr)
    else:
        print(f"SCEPTRE worker completed: {detail}; details: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
