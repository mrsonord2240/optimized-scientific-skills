import argparse
import numpy as np
import pandas as pd
from pathlib import Path

def median_center(log_int):
    '''Subtract each sample median log2 intensity (corrects LOCATION only, cannot manufacture variance).'''
    sample_medians = log_int.median(axis=0)
    return log_int - sample_medians + sample_medians.median()


def sample_loading_normalize(plex):
    '''Per-channel scalar correcting total load WITHIN one plex; target = mean column sum.'''
    target = plex.sum(axis=0).mean()
    return plex * (target / plex.sum(axis=0))


def irs_scale(plexes, ref_cols):
    '''Internal Reference Scaling bridge (Plubell 2017): pin each plex reference channel to a common per-protein value.'''
    refs = pd.concat([p[ref] for p, ref in zip(plexes, ref_cols)], axis=1)
    refs = refs.where(refs > 0)    # a 0 or missing reference cannot anchor the bridge
    unbridged = refs.index[refs.isna().any(axis=1)]
    if len(unbridged):
        print(f'IRS: {len(unbridged)} proteins lack a reference in >=1 plex; set to NaN: {list(unbridged[:10])}')
    geomean = np.exp(np.log(refs).mean(axis=1, skipna=False))    # per-protein geometric mean across plexes
    out = []
    for i, p in enumerate(plexes):
        factor = geomean / refs.iloc[:, i]    # per-protein per-plex scaling factor
        out.append(p.mul(factor, axis=0))
    return out


def plex_offset(a, b, ref_a, ref_b):
    '''Median per-protein log2 offset between plexes on NON-reference channels (IRS equalizes references by construction).'''
    la = np.log2(a.drop(columns=ref_a)).mean(axis=1)
    lb = np.log2(b.drop(columns=ref_b)).mean(axis=1)
    return float(np.nanmedian(lb - la))


parser = argparse.ArgumentParser()
sub = parser.add_subparsers(dest='mode', required=True)
median = sub.add_parser('median')
median.add_argument('input')
median.add_argument('output')
irs = sub.add_parser('irs')
irs.add_argument('--plex', action='append', required=True)
irs.add_argument('--ref', action='append', required=True)
irs.add_argument('--out-dir', required=True)
args = parser.parse_args()
if args.mode == 'median':
    matrix = pd.read_csv(args.input, index_col=0)
    median_center(np.log2(matrix.replace(0, np.nan))).to_csv(args.output)
else:
    if len(args.plex) != len(args.ref):
        parser.error('provide one --ref per --plex')
    plexes = [pd.read_csv(path, index_col=0) for path in args.plex]
    bridged = irs_scale([sample_loading_normalize(p) for p in plexes], args.ref)
    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for i, frame in enumerate(bridged, 1):
        frame.to_csv(output_dir / f'plex_{i}.csv')
