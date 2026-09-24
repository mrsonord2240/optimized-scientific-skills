"""Aggregate MAGeCK sgRNA-level LFC to per-variant scores, separating bystander-free sgRNAs.

Purpose: per-variant fitness for a base-editing screen after efficiency filtering.
Inputs:  MAGeCK sgrna_summary.txt (lowercase 'sgrna' column) and a per-sgRNA annotation TSV with columns
         sgrna, target_variant, n_bystanders.
Output:  <prefix>.target_only.tsv (mean/std/count of LFC per variant, sgRNAs with 0 bystanders) and
         <prefix>.mixed.tsv (sgRNAs with bystanders, for downstream attribution).
Usage:   python aggregate_variant_scores.py sgrna_summary.txt variant_annotation.tsv --out-prefix out/variant_scores
Needs:   pandas.
"""
import pandas as pd


def aggregate_variant_scores(mageck_sgrna_summary, variant_annotation_df):
    '''Aggregate sgRNA-level scores to per-variant scores.
    variant_annotation_df: per-sgRNA -> predicted variants (target + bystanders),
    keyed on the same sgRNA-identifier column name as mageck_sgrna_summary.
    MAGeCK's real sgrna_summary.txt column is lowercase 'sgrna' (not 'sgRNA') --
    build variant_annotation_df with that same column name.'''
    for _name, _frame in (('mageck_sgrna_summary', mageck_sgrna_summary), ('variant_annotation_df', variant_annotation_df)):
        if 'sgrna' not in _frame.columns:
            raise ValueError(f"{_name} is missing the 'sgrna' merge column (got {list(_frame.columns)})")
    df = mageck_sgrna_summary.merge(variant_annotation_df, on='sgrna')
    # Target-only contribution: sgRNAs with no bystanders
    target_only = df[df['n_bystanders'] == 0]
    target_only_scores = target_only.groupby('target_variant')['LFC'].agg(['mean', 'std', 'count'])
    # Mixed signal: sgRNAs with bystanders
    mixed = df[df['n_bystanders'] > 0]
    return target_only_scores, mixed


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('sgrna_summary', help='MAGeCK *.sgrna_summary.txt')
    ap.add_argument('variant_annotation', help='TSV with columns sgrna, target_variant, n_bystanders')
    ap.add_argument('--out-prefix', default='variant_scores', help='writes <prefix>.target_only.tsv and <prefix>.mixed.tsv')
    a = ap.parse_args()
    summ = pd.read_csv(a.sgrna_summary, sep='\t')
    ann = pd.read_csv(a.variant_annotation, sep='\t')
    target_only_scores, mixed = aggregate_variant_scores(summ, ann)
    target_only_scores.to_csv(a.out_prefix + '.target_only.tsv', sep='\t')
    mixed.to_csv(a.out_prefix + '.mixed.tsv', sep='\t', index=False)
    print(target_only_scores.to_string())
    print(f'{len(mixed)} sgRNAs with bystanders written to {a.out_prefix}.mixed.tsv')
