#!/usr/bin/env python3
'''Decision-curve net benefit of predicted probabilities against treat-all and treat-none.

Net benefit = TP/n - (FP/n) * pt/(1-pt), at threshold probability pt.

Inputs:  CSV with a binary label column and a calibrated predicted-probability column from held-out samples.
Usage:   python decision_curve.py predictions.csv --label label --prob prob [--thresholds 0.05,0.1,0.2,0.3]
Output:  one row per threshold: model net benefit, treat-all net benefit, treat-none (0), and whether the
         model beats both.
'''
import argparse

import numpy as np
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('csv')
    parser.add_argument('--label', required=True)
    parser.add_argument('--prob', required=True)
    parser.add_argument('--thresholds', default='0.05,0.1,0.2,0.3,0.4,0.5')
    args = parser.parse_args()

    data = pd.read_csv(args.csv)
    y, p = data[args.label].to_numpy(), data[args.prob].to_numpy(float)
    if set(np.unique(y)) - {0, 1}:
        raise SystemExit('label must be 0/1')
    if p.min() < 0 or p.max() > 1:
        raise SystemExit('prob must be probabilities in [0, 1], not scores or logits')
    thresholds = [float(t) for t in args.thresholds.split(',')]
    if any(not 0 < t < 1 for t in thresholds):
        raise SystemExit('thresholds must lie strictly between 0 and 1')

    n, prevalence = len(y), y.mean()
    print(f'{args.csv}: {n} rows, prevalence {prevalence:.3f}')
    print('threshold  model   treat-all  treat-none  model beats both')
    for pt in thresholds:
        odds = pt / (1 - pt)
        treat = p >= pt
        model = ((treat & (y == 1)).sum() - (treat & (y == 0)).sum() * odds) / n
        treat_all = prevalence - (1 - prevalence) * odds
        print(f'{pt:9.2f}  {model:6.3f}  {treat_all:9.3f}  {0:10.3f}  {"yes" if model > max(treat_all, 0) else "no"}')
