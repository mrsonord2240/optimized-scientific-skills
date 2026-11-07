"""Aalen-Johansen cumulative incidence from a research outcome table."""
import argparse
import numpy as np
import pandas as pd
from lifelines import AalenJohansenFitter

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("time_column")
p.add_argument("cause_column")
p.add_argument("event_of_interest", type=int)
p.add_argument("output_csv")
args = p.parse_args()
d = pd.read_csv(args.input)
if not {args.time_column, args.cause_column}.issubset(d.columns):
    p.error("time or cause column missing")
time = pd.to_numeric(d[args.time_column], errors="raise")
cause = pd.to_numeric(d[args.cause_column], errors="raise")
if not np.isfinite(time).all() or (time < 0).any() or not np.isfinite(cause).all() or (cause < 0).any() or (cause % 1 != 0).any() or args.event_of_interest <= 0 or not (cause == args.event_of_interest).any():
    p.error("invalid follow-up or integer event code")
fit = AalenJohansenFitter(seed=0).fit(time, cause.astype(int), event_of_interest=args.event_of_interest)
fit.cumulative_density_.to_csv(args.output_csv)
