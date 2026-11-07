import argparse
import numpy as np
import pandas as pd
from scipy import stats

R2_MIN = 0.95
TRANSITION_KEY = ['ModifiedSequence', 'PrecursorCharge', 'FragmentType', 'FragmentSeriesNumber', 'FragmentCharge']

def calibrate_irt(anchor_irt, observed_rt):
    slope, intercept, r, _, _ = stats.linregress(anchor_irt, observed_rt)
    if r ** 2 < R2_MIN:
        raise ValueError(f'iRT fit R^2={r**2:.3f} < {R2_MIN}; gradient may be nonlinear, use LOWESS')
    return slope, intercept, r ** 2

def merge_libraries(libs):
    combined = pd.concat(libs, ignore_index=True)
    combined['precursor_total'] = combined.groupby(['ModifiedSequence', 'PrecursorCharge'])['LibraryIntensity'].transform('sum')
    combined = combined.sort_values('precursor_total', ascending=False)
    return combined.drop_duplicates(subset=TRANSITION_KEY).drop(columns='precursor_total')

def library_stats(lib):
    n_prec = lib.groupby(['ModifiedSequence', 'PrecursorCharge']).ngroups
    return {'precursors': n_prec, 'proteins': lib['ProteinId'].nunique(),
            'transitions_per_precursor': round(len(lib) / n_prec, 1)}

parser = argparse.ArgumentParser()
sub = parser.add_subparsers(dest='mode', required=True)
cal = sub.add_parser('calibrate')
cal.add_argument('anchors')
cal.add_argument('library')
cal.add_argument('output')
merge = sub.add_parser('merge')
merge.add_argument('output')
merge.add_argument('libraries', nargs='+')
args = parser.parse_args()
if args.mode == 'calibrate':
    anchors = pd.read_csv(args.anchors, sep='\t')
    library = pd.read_csv(args.library, sep='\t')
    slope, intercept, r2 = calibrate_irt(anchors['irt'], anchors['observed_rt'])
    library['calibrated_rt'] = slope * library['iRT'] + intercept
    library.to_csv(args.output, sep='\t', index=False)
    print(f'slope={slope}, intercept={intercept}, r2={r2}')
else:
    merged = merge_libraries([pd.read_csv(path, sep='\t') for path in args.libraries])
    merged.to_csv(args.output, sep='\t', index=False)
    print(library_stats(merged))
