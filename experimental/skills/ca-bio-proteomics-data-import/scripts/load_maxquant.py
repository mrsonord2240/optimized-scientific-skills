import pandas as pd
import numpy as np
def load_clean_maxquant(path):
    pg = pd.read_csv(path, sep='\t', low_memory=False)  # mixed-type columns
    if 'Protein IDs' not in pg:
        raise ValueError("Expected MaxQuant proteinGroups.txt with a 'Protein IDs' column; this is not that table")
    flag_cols = ('Reverse', 'Potential contaminant', 'Only identified by site')
    if not any(c in pg for c in flag_cols):
        raise ValueError("Expected MaxQuant proteinGroups.txt bookkeeping columns; received none of Reverse, Potential contaminant, Only identified by site")
    keep = pd.Series(True, index=pg.index)
    for col in flag_cols:
        if col in pg:
            keep &= pg[col].fillna('').ne('+')
    rows_read, rows_kept = len(pg), int(keep.sum())
    pg = pg.loc[keep].copy()

    # Semicolon lists: first entry is the leading/razor identifier; Gene names may be blank
    pg['leading_protein'] = pg['Protein IDs'].str.split(';').str[0]
    pg['leading_gene'] = pg.get('Gene names', pd.Series('', index=pg.index)).fillna('').str.split(';').str[0]

    lfq_cols = [c for c in pg.columns if c.startswith('LFQ intensity ')]  # MaxLFQ-normalized, between-sample comparable
    if not lfq_cols:
        tmt_cols = [c for c in pg.columns if c.startswith('Reporter intensity corrected ')]
        if tmt_cols:
            raise ValueError("No LFQ intensity columns: this looks like MaxQuant TMT. Use 'Reporter intensity corrected <channel>' columns; reporter processing -> quantification")
        raise ValueError('No LFQ intensity columns: use a label-free MaxQuant proteinGroups.txt, or verify the supported quant family')
    matrix = pg[['leading_protein', 'leading_gene'] + lfq_cols].copy()
    matrix[lfq_cols] = matrix[lfq_cols].replace(0, np.nan)  # 0 means not-quantified; log2(0) = -inf
    matrix[lfq_cols] = np.log2(matrix[lfq_cols])
    matrix = matrix[matrix[lfq_cols].notna().any(axis=1)]  # drop groups with no valid LFQ value
    return matrix, lfq_cols, {'rows_read': rows_read, 'after_bookkeeping': rows_kept, 'quantified': len(matrix)}

def assess_missingness(matrix, sample_cols):
    total_pct = 100 * matrix[sample_cols].isna().sum().sum() / matrix[sample_cols].size
    mean_abund = matrix[sample_cols].mean(axis=1)  # negative corr with missingness => MNAR / left-censored
    mnar_corr = mean_abund.corr(matrix[sample_cols].isna().sum(axis=1))
    return total_pct, mnar_corr

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    args = parser.parse_args()
    matrix, columns, counts = load_clean_maxquant(args.input)
    matrix.to_csv(args.output, sep='\t', index=False)
    print(counts)
    print('missingness percentage and abundance correlation:', assess_missingness(matrix, columns))
