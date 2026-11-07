"""Partition CRISPResso2 alleles by edit pattern at the target base and its bystanders.

Purpose: bystander attribution -- how much of the read population carries target-only,
target+bystander, or bystander-only edits.
Inputs:  Alleles_frequency_table.zip, 1-indexed target position and bystander positions
         (positions within the aligned amplicon string).
Output:  table of target_edited x bystander_<pos>_edited with summed %Reads.
Usage:   python deconvolute_bystander.py Alleles_frequency_table.zip --target-pos 67 --bystander-pos 69
Needs:   pandas.
"""
import pandas as pd


def deconvolute_bystander(allele_table_path, target_pos, bystander_pos_list):
    '''From CRISPResso2 allele table, partition reads by edit pattern at target + bystanders.
    Returns: per-pattern frequency for each combination of target/bystander edits.'''
    alleles = pd.read_csv(allele_table_path, sep='\t', compression='zip')
    required_cols = {'Aligned_Sequence', 'Reference_Sequence', '%Reads'}
    missing = required_cols - set(alleles.columns)
    if missing:
        raise ValueError(f"Unexpected Alleles_frequency_table schema: missing {missing}; "
                          f"got columns {list(alleles.columns)}")
    # Mark target_edited and per-bystander_edited
    alleles['target_edited'] = alleles['Aligned_Sequence'].str[target_pos-1] != alleles['Reference_Sequence'].str[target_pos-1]
    for bp in bystander_pos_list:
        alleles[f'bystander_{bp}_edited'] = alleles['Aligned_Sequence'].str[bp-1] != alleles['Reference_Sequence'].str[bp-1]
    # Real Alleles_frequency_table.zip has no 'Reference_pct' column -- the
    # per-allele read-fraction column is '%Reads' (verified against actual
    # CRISPResso2 2.3.4 output).
    return alleles.groupby(['target_edited'] + [f'bystander_{bp}_edited' for bp in bystander_pos_list])['%Reads'].sum().reset_index()


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('allele_table', help='CRISPResso2 Alleles_frequency_table.zip')
    ap.add_argument('--target-pos', type=int, required=True, help='1-indexed position in the aligned allele string')
    ap.add_argument('--bystander-pos', type=int, nargs='*', default=[], help='1-indexed bystander positions')
    ap.add_argument('--out', help='write TSV here (default: print)')
    a = ap.parse_args()
    res = deconvolute_bystander(a.allele_table, a.target_pos, a.bystander_pos)
    if a.out:
        res.to_csv(a.out, sep='\t', index=False)
    else:
        print(res.to_string(index=False))
