"""AP-MS scoring against NEGATIVE-CONTROL IPs (not the input lysate); no median/SL/IRS normalization.

Input : CSV, prey x replicate RAW intensities or spectral counts (first column = prey id); 0/NaN = not detected.
Output: per-prey n_bait, n_ctrl, n_ctrl_runs_used, log2_enrichment, worst_case_log2, interactor (CSV with --out).
Usage : python scripts/apms_score.py ip.csv --bait Bait1,Bait2,Bait3 --ctrl Ctrl1,Ctrl2,Ctrl3 [--min-bait-reps 2] [--fc-cutoff 2.0] [--out scores.csv]
Import: from apms_score import score_vs_control_ips
Checked: numpy 2.5.3, pandas 3.0.5
"""
import numpy as np, pandas as pd

def score_vs_control_ips(ip, bait_cols, ctrl_cols, fc_cutoff=2.0, min_bait_reps=None):
    '''ip: prey x replicate RAW intensities or spectral counts; 0/NaN = not detected. Prey never
    seen in any bait IP get a NaN enrichment (nothing was measured) and are never called.'''
    L = np.log2(ip[list(bait_cols) + list(ctrl_cols)].replace(0, np.nan).astype(float))
    if min_bait_reps is None:
        min_bait_reps = len(bait_cols)          # default: every bait replicate, reproducibility first
    n_bait, n_ctrl = L[bait_cols].notna().sum(axis=1), L[ctrl_cols].notna().sum(axis=1)
    # a control IP that produced no data is skipped by every mean/min below, so say so out loud
    dead = [c for c in ctrl_cols if L[c].notna().sum() == 0]
    if dead:
        print(f'AP-MS: control runs with no data, excluded: {dead}')
    if len(dead) == len(ctrl_cols):
        raise ValueError('every control IP is empty -- nothing to score against')
    Lf = L.copy()
    Lf[ctrl_cols] = Lf[ctrl_cols].fillna(L[ctrl_cols].min())      # per-control-run detection floor
    enrich = Lf[bait_cols].mean(axis=1) - Lf[ctrl_cols].mean(axis=1)
    worst = Lf[bait_cols].min(axis=1) - Lf[ctrl_cols].max(axis=1)  # weakest bait rep vs best control
    out = pd.DataFrame({'n_bait': n_bait, 'n_ctrl': n_ctrl, 'n_ctrl_runs_used': len(ctrl_cols) - len(dead),
                        'log2_enrichment': enrich, 'worst_case_log2': worst})
    out['interactor'] = (n_bait >= min_bait_reps) & ((n_ctrl == 0) | (enrich >= fc_cutoff))
    return out.sort_values('log2_enrichment', ascending=False)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('ip_csv')
    ap.add_argument('--bait', required=True, help='comma-separated bait replicate columns')
    ap.add_argument('--ctrl', required=True, help='comma-separated control IP columns')
    ap.add_argument('--min-bait-reps', type=int, default=None)
    ap.add_argument('--fc-cutoff', type=float, default=2.0)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    scores = score_vs_control_ips(pd.read_csv(a.ip_csv, index_col=0), a.bait.split(','), a.ctrl.split(','),
                                  fc_cutoff=a.fc_cutoff, min_bait_reps=a.min_bait_reps)
    print(f"called {int(scores['interactor'].sum())} of {len(scores)} prey")
    if a.out:
        scores.to_csv(a.out)
