"""Recount preserved results without models, Docker, or third-party packages."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent


def counts(rows):
    completed = [r for r in rows if not r.get("infra_error", False)]
    return {
        "attempts": len(rows),
        "scored": len(completed),
        "passed": sum(bool(r["passed"]) for r in completed),
        "any_skill": sum(bool(r.get("skills_touched")) for r in completed),
        "infra_errors": len(rows) - len(completed),
    }


def summarize():
    records = [json.loads(line) for line in (HERE / "results.jsonl").read_text(
        encoding="utf-8").splitlines() if line.strip()]
    groups = defaultdict(list)
    sources = set()
    for item in records:
        if item["source"] in sources:
            raise ValueError(f"duplicate result source: {item['source']}")
        sources.add(item["source"])
        r = item["result"]
        runset = Path(item["source"]).parts[1]
        key = (runset, r["model"], r["arm"], r.get("bundle", "all"), r.get("desc", "full"))
        groups[key].append(r)
    experiments = []
    for (runset, model, arm, bundle, desc), rows in sorted(groups.items()):
        experiments.append({
            "runset": runset, "model": model, "arm": arm,
            "bundle": bundle, "description": desc, **counts(rows),
            "tasks": {task: counts([r for r in rows if r["task"] == task])
                      for task in sorted({r["task"] for r in rows})},
        })
    matched = [x["result"] for x in records
               if x["source"].startswith("runner/runs_b12_short/")]
    tasks = sorted({r["task"] for r in matched})
    comparison = {}
    for runset in ("runs_base", "runs_b12_full", "runs_b12_short"):
        rows = [x["result"] for x in records
                if x["source"].startswith(f"runner/{runset}/")
                and x["result"]["task"] in tasks]
        cells = {(r["task"], r["rep"]) for r in rows}
        expected = {(task, rep) for task in tasks for rep in range(1, 6)}
        if cells != expected or len(rows) != len(expected):
            raise ValueError(f"missing or duplicate matched repeat in {runset}")
        comparison[runset] = {
            **counts(rows),
            "tasks": {task: counts([r for r in rows if r["task"] == task]) for task in tasks},
        }
    return {"result_count": len(records), "matched_tasks": tasks,
            "selective_comparison": comparison, "experiments": experiments}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="compare with the committed summary")
    args = parser.parse_args()
    summary = summarize()
    if args.check:
        saved = json.loads((HERE / "summary.json").read_text(encoding="utf-8"))
        if summary != saved:
            raise SystemExit("FAIL: summary.json differs from the preserved results")
        print(f"PASS: {summary['result_count']} results; all seven-task comparison cells complete")
    else:
        print(json.dumps(summary, indent=2) + "\n", end="")


if __name__ == "__main__":
    main()
