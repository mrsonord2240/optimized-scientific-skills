# Reference: pandas 2.2+, numpy 1.26+, scipy 1.12+ | Verify API if version differs
#
# Genetic-interaction (GI) scoring for combinatorial CRISPR screens.
# Identifies synthetic-lethal and synthetic-rescue paralog pairs.

import pandas as pd
import numpy as np
from scipy.stats import zscore, norm
from statsmodels.stats.multitest import multipletests

# === INPUTS ===
# paired_lfc.tsv: cassette_id, gene_A, gene_B, lfc (paired double-KO LFC vs control)
# single_lfc.tsv: gene, lfc (single-KO LFC vs control)
paired_df = pd.read_csv('paired_lfc.tsv', sep='\t')
single_df = pd.read_csv('single_lfc.tsv', sep='\t')

# === STEP 1: BUILD SINGLE-GENE LOOKUP ===
single_lookup = dict(zip(single_df['gene'], single_df['lfc']))

# === STEP 2: AGGREGATE CASSETTE-LEVEL TO PAIR-LEVEL ===
# Inzolia has 4 cassettes per pair; aggregate mean LFC across replicates
pair_lfc = paired_df.groupby(['gene_A', 'gene_B']).agg(
    paired_lfc_mean=('lfc', 'mean'),
    paired_lfc_std=('lfc', 'std'),
    n_cassettes=('lfc', 'count')
).reset_index()

# === STEP 3: COMPUTE EXPECTED ADDITIVE FROM SINGLETONS ===
pair_lfc['single_A_lfc'] = pair_lfc['gene_A'].map(single_lookup)
pair_lfc['single_B_lfc'] = pair_lfc['gene_B'].map(single_lookup)
pair_lfc['expected_additive'] = pair_lfc['single_A_lfc'] + pair_lfc['single_B_lfc']

# Drop pairs where singleton LFC is missing
pair_lfc = pair_lfc.dropna(subset=['single_A_lfc', 'single_B_lfc'])

# === STEP 4: COMPUTE GI SCORE ===
# GI = observed_paired - expected_additive
# Negative GI = synthetic lethal (paired more depleted than expected)
# Positive GI = synthetic rescue (paired less depleted than expected)
pair_lfc['gi_score'] = pair_lfc['paired_lfc_mean'] - pair_lfc['expected_additive']

# === STEP 5: Z-NORMALIZE GI ===
pair_lfc['gi_z'] = zscore(pair_lfc['gi_score'])

# === STEP 6: CLASSIFY (raw cutoff) ===
# Adequate only for a small, hand-curated pair set (see SKILL.md's Minimum
# pair count note). At genome scale it is NOT FDR-controlled: it only looks
# safe when a few very large-effect true hits inflate the population SD
# enough to suppress noise from crossing the threshold. See STEP 6.5.
pair_lfc['gi_class'] = np.where(pair_lfc['gi_z'] < -2, 'synthetic_lethal',
                                 np.where(pair_lfc['gi_z'] > 2, 'synthetic_rescue',
                                          'no_interaction'))

# === STEP 6.5: FDR-CORRECT FOR GENOME SCALE ===
# Convert each z to a two-sided normal-tail p-value and apply Benjamini-Hochberg
# across all tested pairs. Checked on synthetic data (Inzolia scale, 4,435
# pairs): an all-null screen (no true interactions) called ~200 pairs
# "significant" by the raw cutoff alone (chance, ~4.5%) vs. 0 after BH; 20
# true hits at a modest, realistic effect size (GI=-0.3, SNR~5) gave 178 raw
# calls of which 158 (89%) were false, vs. 0 false positives (17/20 true
# hits recovered) after BH. Use this classification for genome-scale calls;
# use STEP 6's raw classification only for the small hand-curated case.
pair_lfc['gi_pvalue'] = 2 * norm.sf(np.abs(pair_lfc['gi_z']))
pair_lfc['gi_fdr_reject'], pair_lfc['gi_fdr'], _, _ = multipletests(
    pair_lfc['gi_pvalue'], alpha=0.05, method='fdr_bh')
pair_lfc['gi_class_fdr'] = np.where(
    pair_lfc['gi_fdr_reject'] & (pair_lfc['gi_z'] < 0), 'synthetic_lethal',
    np.where(pair_lfc['gi_fdr_reject'] & (pair_lfc['gi_z'] > 0), 'synthetic_rescue',
             'no_interaction'))

# === STEP 7: OUTPUT ===
# Synthetic-lethal candidates (drug-target combinations), FDR-corrected
sl_pairs = pair_lfc[pair_lfc['gi_class_fdr'] == 'synthetic_lethal'].sort_values('gi_z')
print(f"Synthetic-lethal pairs, raw z<-2 (uncorrected): {(pair_lfc['gi_class'] == 'synthetic_lethal').sum()}")
print(f'Synthetic-lethal pairs, BH-FDR<0.05 (recommended at genome scale): {len(sl_pairs)}')
print(sl_pairs[['gene_A', 'gene_B', 'paired_lfc_mean', 'expected_additive',
                 'gi_score', 'gi_z', 'gi_fdr']].head(20).to_string(index=False))

# Synthetic-rescue (compensatory pathways), FDR-corrected
sr_pairs = pair_lfc[pair_lfc['gi_class_fdr'] == 'synthetic_rescue'].sort_values('gi_z', ascending=False)
print(f"\nSynthetic-rescue pairs, raw z>2 (uncorrected): {(pair_lfc['gi_class'] == 'synthetic_rescue').sum()}")
print(f'Synthetic-rescue pairs, BH-FDR<0.05 (recommended at genome scale): {len(sr_pairs)}')

# === EXPORT ===
pair_lfc.to_csv('gi_scores.tsv', sep='\t', index=False)
sl_pairs.to_csv('synthetic_lethal_pairs.tsv', sep='\t', index=False)
sr_pairs.to_csv('synthetic_rescue_pairs.tsv', sep='\t', index=False)

# === VALIDATION ===
# For top synthetic-lethal hits, cross-validate against:
# 1. Known paralog pairs (e.g., MAPK1/MAPK3, AKT1/AKT2)
# 2. Multiple cell lines (Dede 2020: 19 of 24 (79%) SL pairs reproduce in at least 2 of 3 lines; 14 of 24 (58%) in all 3)
# 3. Orthogonal modality (CRISPRi if originally Cas9, or vice versa)
# 4. Arrayed validation in target cell line
known_paralog_pairs = [('MAPK1', 'MAPK3'), ('AKT1', 'AKT2'),
                        ('PIK3CA', 'PIK3CB'), ('HSP90AA1', 'HSP90AB1')]
recovered = [(a, b) for a, b in known_paralog_pairs
              if ((pair_lfc['gene_A'] == a) & (pair_lfc['gene_B'] == b) &
                  (pair_lfc['gi_class_fdr'] == 'synthetic_lethal')).any()
              or ((pair_lfc['gene_A'] == b) & (pair_lfc['gene_B'] == a) &
                  (pair_lfc['gi_class_fdr'] == 'synthetic_lethal')).any()]
print(f'\nKnown paralog pairs recovered: {len(recovered)}/{len(known_paralog_pairs)}')
