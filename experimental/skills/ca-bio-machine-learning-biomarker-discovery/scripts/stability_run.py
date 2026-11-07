"""Source subsampling stability selection on independent-sample data."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from stability_selection import stability_selection

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("outcome")
p.add_argument("output_csv")
p.add_argument("summary_json")
args = p.parse_args()
data = pd.read_csv(args.input)
if args.outcome not in data:
    p.error("outcome column missing")
y = data.pop(args.outcome)
if y.isna().any() or y.nunique() != 2:
    p.error("binary outcome is required")
X = data.apply(pd.to_numeric, errors="raise")
if X.isna().any().any() or X.shape[1] < 5:
    p.error("complete numeric features required, at least five")
q = max(1, int(np.sqrt((2 * 0.6 - 1) * X.shape[1])))
freq, stability = stability_selection(X, y, q=q, seed=0)
null_freq, _ = stability_selection(X, np.random.default_rng(1).permutation(y), q=q, seed=0)
pd.DataFrame({"feature": X.columns, "frequency": freq, "null_frequency": null_freq}).to_csv(args.output_csv, index=False)
summary = {"q": q, "threshold": 0.6, "stability": None if not np.isfinite(stability) else float(stability),
           "stable_count": int((freq > 0.6).sum()), "null_count": int((null_freq > 0.6).sum())}
Path(args.summary_json).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
