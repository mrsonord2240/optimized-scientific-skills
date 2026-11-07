"""Fit an exploratory fixed-penalty Cox model, not a clinical predictor."""
import argparse
import numpy as np
import pandas as pd
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("time_column")
p.add_argument("event_column")
p.add_argument("output_csv")
p.add_argument("--alpha", type=float, required=True)
p.add_argument("--l1-ratio", type=float, default=0.9)
args = p.parse_args()
data = pd.read_csv(args.input)
if not {args.time_column, args.event_column}.issubset(data.columns) or args.time_column == args.event_column:
    p.error("distinct time and event columns required")
time = pd.to_numeric(data.pop(args.time_column), errors="raise")
event = pd.to_numeric(data.pop(args.event_column), errors="raise")
X = data.apply(pd.to_numeric, errors="raise")
if X.empty or not np.isfinite(X.to_numpy()).all() or not np.isfinite(time).all() or (time < 0).any() or not event.isin([0, 1]).all() or not (event == 1).any():
    p.error("invalid features, time or 0/1 event")
if not np.isfinite(args.alpha) or args.alpha <= 0 or not 0 < args.l1_ratio <= 1:
    p.error("alpha must be positive and l1-ratio in (0,1]")
y = Surv.from_arrays(event=event.astype(bool), time=time.astype(float))
model = CoxnetSurvivalAnalysis(l1_ratio=args.l1_ratio, alphas=[args.alpha]).fit(X, y)
pd.DataFrame({"feature": X.columns, "coefficient": model.coef_[:, 0]}).to_csv(args.output_csv, index=False)
