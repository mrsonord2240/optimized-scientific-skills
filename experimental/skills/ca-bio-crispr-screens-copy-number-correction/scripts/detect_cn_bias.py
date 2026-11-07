import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

def detect_cn_bias(gene_lfc_df, cn_df, amplified_cn=4, diploid_cn=(1.5, 2.5), low_power_n=8):
    '''Test whether amplified genes are depleted relative to diploid ones.

    The genome-wide Spearman rho alone is not enough: a real focal amplicon covers a handful of
    genes out of ~18,000, so rho stays near zero while those genes are severely depleted. On a
    real HAP1 screen with a planted 8-gene 17q12 amplicon at CN 15, rho was -0.034 (not
    "biased") while amplified genes averaged LFC -2.64 against -0.18 for diploid genes. Report
    the stratified gap, and treat either signal as bias.

    n_amplified>=3 is only the floor for the Mann-Whitney test to run at all, not a guarantee it
    has power to detect a real effect. On the real HT-29 FAM84B/MYC/POU5F1B block (CN=8)
    post-CRISPRcleanR, the residual gap was still -0.98 logFC (MYC itself was left uncorrected;
    its two flanking genes were), but n=3 gave p=0.208 -- bias_present reads False even though
    the bias is real and substantial. Bootstrapping that same real 3-gene mixture at larger n
    shows why: P(p<0.01) is only ~0.04 at n=3, ~0.13 at n=8, ~0.30 at n=15 -- the test needs
    dozens of amplified genes at this effect size before it reliably clears its own threshold.
    Below `low_power_n` amplified genes, do not read bias_present=False as "correction
    succeeded" -- check `suspicious_despite_ns` and the raw `amplified_vs_diploid_gap` instead.
    '''
    merged = gene_lfc_df.merge(cn_df, on='gene')
    rho, p = spearmanr(merged['copy_number'], merged['lfc'])
    amplified = merged[merged['copy_number'] > amplified_cn]['lfc']
    diploid = merged[merged['copy_number'].between(*diploid_cn)]['lfc']
    gap, p_gap = float('nan'), float('nan')
    if len(amplified) >= 3 and len(diploid) >= 3:
        gap = amplified.mean() - diploid.mean()
        p_gap = mannwhitneyu(amplified, diploid, alternative='less').pvalue
    low_power = len(amplified) < low_power_n
    return {
        'cn_lfc_rho': rho,
        'p_value': p,
        'n_amplified_genes': len(amplified),
        'amplified_mean_lfc': amplified.mean(),
        'diploid_mean_lfc': diploid.mean(),
        'amplified_vs_diploid_gap': gap,          # negative = amplified genes more depleted
        'p_amplified_more_depleted': p_gap,
        'bias_present': bool((rho < -0.1 and p < 0.01) or (gap < -0.5 and p_gap < 0.01)),
        'low_power_floor': bool(low_power),        # test ran, but underpowered below low_power_n
        'suspicious_despite_ns': bool(low_power and not pd.isna(gap) and gap < -0.5),
    }

if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description='Assess CN versus gene LFC')
    parser.add_argument('gene_lfc')
    parser.add_argument('copy_number')
    parser.add_argument('output')
    args = parser.parse_args()
    result = detect_cn_bias(pd.read_csv(args.gene_lfc, sep='\t'), pd.read_csv(args.copy_number, sep='\t'))
    with open(args.output, 'w', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
