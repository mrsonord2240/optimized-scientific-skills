#!/usr/bin/env python
"""CEGv2 / NEGv1 essentialome recovery (precision-recall AUC) from gene-level LFC.

Purpose: check the screen has detectable essentiality signal (PR-AUC >0.7 passes; <0.5 = no signal).
Inputs:  gene-level LFC table with columns gene, lfc (negative = depleted); CEGv2.txt and NEGv1.txt
         from https://github.com/hart-lab/bagel (a GENE header column, or one symbol per line).
Usage:   python essentialome_recovery.py gene_lfc.tsv CEGv2.txt NEGv1.txt
Import:  from essentialome_recovery import essentialome_recovery
"""
import argparse
import pandas as pd
from sklearn.metrics import precision_recall_curve, auc, roc_auc_score


def essentialome_recovery(gene_lfc_df, cegv2_set, negv1_set):
    '''gene_lfc_df: must have ["gene", "lfc"] columns (gene-level mean LFC, negative = depleted).
    cegv2_set, negv1_set: sets of gene symbols from Hart 2017.'''
    labeled = gene_lfc_df[gene_lfc_df['gene'].isin(cegv2_set | negv1_set)].copy()
    labeled['is_essential'] = labeled['gene'].isin(cegv2_set).astype(int)
    y_score = -labeled['lfc']  # negative LFC = depleted = more essential -> higher score
    precision, recall, _ = precision_recall_curve(labeled['is_essential'], y_score)
    return {
        'pr_auc': auc(recall, precision),
        'roc_auc': roc_auc_score(labeled['is_essential'], y_score),
        'n_essential_detected': labeled['is_essential'].sum(),
        'n_nonessential_detected': (1 - labeled['is_essential']).sum(),
    }


def read_gene_set(path):
    '''BAGEL reference lists have a GENE header (plus other columns); otherwise one symbol per line.'''
    df = pd.read_csv(path, sep='\t')
    if 'GENE' in df.columns:
        return set(df['GENE'])
    return set(pd.read_csv(path, sep='\t', header=None)[0])


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('gene_lfc', help='TSV with columns gene, lfc')
    ap.add_argument('cegv2', help='CEGv2.txt')
    ap.add_argument('negv1', help='NEGv1.txt')
    args = ap.parse_args()
    res = essentialome_recovery(pd.read_csv(args.gene_lfc, sep='\t'),
                                read_gene_set(args.cegv2), read_gene_set(args.negv1))
    for k, v in res.items():
        print(f'{k}: {v:.4f}' if isinstance(v, float) else f'{k}: {v}')
    print('PASS (PR-AUC > 0.7)' if res['pr_auc'] > 0.7 else
          'NO SIGNAL (PR-AUC < 0.5)' if res['pr_auc'] < 0.5 else 'CONCERNING (0.5 <= PR-AUC <= 0.7)')
