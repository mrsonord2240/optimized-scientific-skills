import pandas as pd
import numpy as np

def time_course_consistency(mle_results, conditions=['day7', 'day14', 'day21']):
    '''Identify genes with monotonic beta trends across timepoints.
    Returns genes where all betas same sign and trend is monotone.'''
    beta_cols = [f'{c}|beta' for c in conditions]
    fdr_cols = [f'{c}|fdr' for c in conditions]
    df = mle_results[['Gene'] + beta_cols + fdr_cols].copy()
    df['all_negative'] = (df[beta_cols] < 0).all(axis=1)
    df['all_positive'] = (df[beta_cols] > 0).all(axis=1)
    df['monotone'] = df[beta_cols].apply(lambda x: (np.diff(x) <= 0).all() or (np.diff(x) >= 0).all(), axis=1)
    df['any_sig'] = (df[fdr_cols] < 0.05).any(axis=1)
    return df[df['monotone'] & df['any_sig']].sort_values(beta_cols[-1])

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--conditions', default='day7,day14,day21')
    args = parser.parse_args()
    time_course_consistency(pd.read_csv(args.input, sep='\t'), args.conditions.split(',')).to_csv(args.output, sep='\t', index=False)
