"""Summarise <runs dir>/*/result.json: pass rate, uptake and cost per task x condition x model.

  python report.py runs_base runs_b12_full runs_b12_short    # several run dirs in one table

A condition is the arm plus, for skill arms, the bundle and description mode the run recorded.

Infra errors (API failures) are excluded from pass rates and counted separately, so a
network blip is never scored as the model failing.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = [Path(a) for a in sys.argv[1:]] or [HERE / "runs"]

sys.path.insert(0, str(HERE))
from run import TASKS, relevant_skills  # noqa: E402

RELEVANT = {t: set(relevant_skills(t)) for t in TASKS}

groups = defaultdict(list)
totals = defaultdict(list)
for f in sorted(f for d in RUNS for f in d.glob("*/result.json")):
    r = json.loads(f.read_text(encoding="utf-8"))
    cond = r["arm"] if r["arm"] == "none" else f"{r['arm']}/{r.get('bundle', 'all')}/{r.get('desc', 'full')}"
    groups[(r["model"], r["task"], cond)].append(r)

print(f"{'model':32s} {'task':22s} {'condition':22s} {'n':>3s} {'pass':>9s} {'any-skill':>9s} {'rel-skill':>9s} {'$/run':>7s} {'steps':>5s} infra")
for (model, task, arm), rs in sorted(groups.items()):
    ok = [r for r in rs if not r["infra_error"]]
    n = len(ok)
    if not n:
        print(f"{model:32s} {task:22s} {arm:22s}   0 (all {len(rs)} runs infra errors)")
        continue
    p = sum(r["passed"] for r in ok)
    anys = sum(bool(r["skills_touched"]) for r in ok)
    rel = sum(bool(set(r["skills_touched"]) & RELEVANT.get(task, set())) for r in ok)
    cost = sum(r["cost_usd"] for r in ok) / n
    steps = sum(r["steps"] for r in ok) / n
    up = "-" if arm == "none" else f"{anys}/{n}"
    rl = "-" if arm == "none" else f"{rel}/{n}"
    totals[(model, arm)].append((p / n, anys / n, rel / n, cost))
    print(f"{model:32s} {task:22s} {arm:22s} {n:3d} {p:3d}/{n:<3d}{100 * p / n:4.0f}% {up:>9s} {rl:>9s} {cost:7.3f} {steps:5.1f} {len(rs) - n}")

print("\nmean over tasks (each task weighted equally)")
for (model, arm), rows in sorted(totals.items()):
    k = len(rows)
    p, a, rl, c = (sum(col) / k for col in zip(*rows))
    up = "" if arm == "none" else f" any-skill {100 * a:3.0f}% rel-skill {100 * rl:3.0f}%"
    print(f"{model:32s} {arm:22s} tasks {k:2d} pass {100 * p:3.0f}%{up} $/run {c:.3f}")
