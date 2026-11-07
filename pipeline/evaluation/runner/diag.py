"""Diagnostics view over runs/*/result.json.

  python diag.py                 # latency / retry / stall table per model x task x arm
  python diag.py --failures      # every non-pass with its failure class and the agent's own note
  python diag.py --live          # runs in flight and the phase they are stuck in
  python diag.py --runs DIR      # use another runs directory

Answers "was that a model failure or an infrastructure failure, and where did the time go".
"""
import argparse
import json
import statistics as st
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("--runs", default=str(HERE / "runs"))
ap.add_argument("--failures", action="store_true")
ap.add_argument("--live", action="store_true")
ap.add_argument("--stall", type=float, default=120, help="seconds in one phase before --live flags it")
a = ap.parse_args()
RUNS = Path(a.runs)

if a.live:
    now = time.time()
    n = 0
    for d in sorted(RUNS.iterdir()):
        sf = d / "status.json"
        if not d.is_dir() or (d / "result.json").exists() or not sf.exists():
            continue
        s = json.loads(sf.read_text(encoding="utf-8"))
        age = now - s["t"]
        n += 1
        flag = "  <-- STALLED?" if age > a.stall and s["phase"] in ("model_call", "exec", "self_report") else ""
        print(f"{d.name[:80]:80s} step={s['step']:2d} phase={s['phase']:15s} {age:6.0f}s{flag}")
    print(f"{n} run(s) in flight")
    raise SystemExit

runs = [json.loads(f.read_text(encoding="utf-8")) for f in sorted(RUNS.glob("*/result.json"))]

if a.failures:
    for r in runs:
        if r["passed"]:
            continue
        d = r.get("diagnostics") or {}
        print(f"\n### {r['task']} / {r['arm']} / r{r['rep']}  [{r.get('failure_class')}]  steps={r['steps']} {r['seconds']:.0f}s ${r['cost_usd']:.3f}")
        if r.get("answer") is not None:
            print("answer:", json.dumps(r["answer"])[:200])
        flags = []
        if d.get("stalled_calls_gap_over_60s"):
            flags.append(f"model calls silent >60s at steps {d['stalled_calls_gap_over_60s']}")
        if d.get("n_retries"):
            flags.append(f"{d['n_retries']} API retries {d['retry_kinds']}")
        if d.get("cmd_timeouts"):
            flags.append(f"command timeouts at steps {d['cmd_timeouts']}")
        if d.get("truncated_by_length"):
            flags.append(f"replies cut off by length at {d['truncated_by_length']}")
        if d.get("empty_replies"):
            flags.append(f"empty replies at {d['empty_replies']}")
        if d.get("failed_call"):
            flags.append(f"API failure: {[x['kind'] for x in d['failed_call']]}")
        if flags:
            print("infra:", "; ".join(flags))
        note = (r.get("self_report") or "").strip().replace("\n", " | ")
        print("agent's own note:", note[:600] if note else "(none)")
    raise SystemExit

groups = defaultdict(list)
for r in runs:
    groups[(r["model"], r["task"], r["arm"])].append(r)
print(f"{'task':24s} {'arm':9s} {'runs':>4s} {'lat p50':>8s} {'p95':>6s} {'max':>6s} {'model%':>6s} {'exec%':>6s} "
      f"{'retry':>5s} {'stall':>5s} {'cmdTO':>5s} {'empty':>5s} {'len':>4s} {'ctr_s':>5s} {'rsn_tok':>8s}  providers")
for (model, task, arm), rs in sorted(groups.items(), key=lambda kv: kv[0][1:]):
    ds = [r["diagnostics"] for r in rs if r.get("diagnostics")]
    if not ds:
        continue
    lat = [d["latency_p50"] for d in ds if d["latency_p50"] is not None]
    tot = sum(r["seconds"] for r in rs) or 1
    mt = sum(d["model_s_total"] for d in ds) / tot * 100
    et = sum(d["exec_s_total"] for d in ds) / tot * 100
    prov = Counter(p for d in ds for p in d["providers"])
    print(f"{task:24s} {arm:9s} {len(rs):4d} {st.median(lat) if lat else 0:8.1f} "
          f"{max((d['latency_p95'] or 0) for d in ds):6.1f} {max((d['latency_max'] or 0) for d in ds):6.1f} "
          f"{mt:6.0f} {et:6.0f} {sum(d['n_retries'] for d in ds):5d} "
          f"{sum(len(d['stalled_calls_gap_over_60s']) for d in ds):5d} {sum(len(d['cmd_timeouts']) for d in ds):5d} "
          f"{sum(len(d['empty_replies']) for d in ds):5d} {sum(len(d['truncated_by_length']) for d in ds):4d} "
          f"{st.mean(d['container_start_s'] for d in ds):5.0f} {sum(d['reasoning_tokens'] for d in ds):8d}  "
          f"{', '.join(f'{p}:{c}' for p, c in prov.most_common(3))}")

fc = Counter(r.get("failure_class") or "pass" for r in runs)
print("\noutcomes:", dict(fc))
