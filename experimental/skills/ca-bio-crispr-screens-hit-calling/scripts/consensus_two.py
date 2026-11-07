'''Consensus hit calling from multiple methods'''
# Reference: mageck 0.5+, matplotlib 3.8+, numpy 1.26+, pandas 2.2+, scipy 1.12+, statsmodels 0.14+ | Verify API if version differs
import pandas as pd
from scipy.stats import hypergeom

# Load results from different methods
mageck = pd.read_csv('mageck.gene_summary.txt', sep='\t')
bagel = pd.read_csv('bagel_bf.txt', sep='\t')

# Standardize column names
mageck = mageck[['id', 'neg|score', 'neg|fdr']].rename(columns={'id': 'gene'})
bagel = bagel[['GENE', 'BF']].rename(columns={'GENE': 'gene'})

# Merge
merged = mageck.merge(bagel, on='gene', how='outer')

# Call hits per method -- thresholds match SKILL.md's Quantitative Thresholds table
# (FDR<0.05, BF>6). Keep this file in sync with that table, not with any other default.
merged['mageck_hit'] = merged['neg|fdr'] < 0.05
merged['bagel_hit'] = merged['BF'] > 6

# Sanity-check that the two hit sets actually overlap more than chance would predict --
# a hit-set pair with no enrichment is the statistical signature of merging results
# from non-comparable experimental designs (e.g. two different comparisons), not a
# real method disagreement. See references/failure-modes.md: "Consensus across 3 methods
# is empty".
n = len(merged)
a = merged['mageck_hit'].fillna(False)
b = merged['bagel_hit'].fillna(False)
k, K, N = int((a & b).sum()), int(a.sum()), int(b.sum())
if K > 0 and N > 0:
    p_overlap = hypergeom.sf(k - 1, n, K, N)
    if p_overlap > 0.05:
        print(f'WARNING: mageck_hit vs bagel_hit overlap not enriched above chance '
              f'(observed={k}, expected~{K * N / n:.1f}, p={p_overlap:.3f}) -- check both '
              f'files came from the SAME experimental comparison before trusting consensus.')

# Consensus
merged['n_methods'] = merged['mageck_hit'].fillna(False).astype(int) + merged['bagel_hit'].fillna(False).astype(int)

# Results
print('=== Consensus Hit Calling ===')
print(f'MAGeCK hits (FDR<0.05): {merged["mageck_hit"].sum()}')
print(f'BAGEL2 hits (BF>6): {merged["bagel_hit"].sum()}')
print(f'Consensus hits (both): {(merged["n_methods"] == 2).sum()}')

# High confidence hits
high_conf = merged[merged['n_methods'] == 2].sort_values('neg|score')
print('\nTop consensus hits:')
print(high_conf[['gene', 'neg|score', 'neg|fdr', 'BF']].head(20).to_string(index=False))

# Save
high_conf.to_csv('consensus_hits.csv', index=False)
