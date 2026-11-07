#!/usr/bin/env python3
'''Calibration of predicted probabilities on held-out samples.

Inputs:  CSV with a binary label column and a predicted-probability column, from samples the model was not
         trained or recalibrated on.
Usage:   python calibration.py predictions.csv --label label --prob prob [--bins 10]
Output:  reliability curve (equal-mass bins), Brier score vs the prevalence-only Brier, and the calibration
         intercept and slope from a logistic fit of the label on logit(prob).
'''
import argparse

import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('csv')
    parser.add_argument('--label', required=True)
    parser.add_argument('--prob', required=True)
    parser.add_argument('--bins', type=int, default=10)
    args = parser.parse_args()

    data = pd.read_csv(args.csv)
    y, p = data[args.label].to_numpy(), data[args.prob].to_numpy(float)
    if set(np.unique(y)) - {0, 1} or len(set(y)) < 2:
        raise SystemExit('label must be 0/1 with both classes present')
    if p.min() < 0 or p.max() > 1:
        raise SystemExit('prob must be probabilities in [0, 1], not scores or logits')

    prevalence = y.mean()
    print(f'{args.csv}: {len(y)} rows, prevalence {prevalence:.3f}, mean predicted {p.mean():.3f}')
    observed, predicted = calibration_curve(y, p, n_bins=args.bins, strategy='quantile')
    print('Reliability (mean predicted -> observed fraction):')
    for pred, obs in zip(predicted, observed):
        print(f'  {pred:.3f} -> {obs:.3f}')
    print(f'Brier = {brier_score_loss(y, p):.3f} (prevalence-only reference {prevalence * (1 - prevalence):.3f})')

    eps = 1e-6
    logit = np.log(np.clip(p, eps, 1 - eps) / (1 - np.clip(p, eps, 1 - eps))).reshape(-1, 1)
    fit = LogisticRegression(C=1e6, max_iter=5000).fit(logit, y)
    print(f'Calibration slope = {fit.coef_[0][0]:.2f} (1 is ideal; below 1 is overconfident), '
          f'intercept = {fit.intercept_[0]:.2f} (0 is ideal)')
