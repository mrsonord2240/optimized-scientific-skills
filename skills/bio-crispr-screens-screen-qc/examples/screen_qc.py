'''CRISPR screen quality control'''
# Reference: matplotlib 3.8+, numpy 1.26+, pandas 2.2+, seaborn 0.13+ | Verify API if version differs
# Thresholds and pairing follow SKILL.md: stage-specific Gini/zero-count bands, and replicate
# correlation only between samples of the same condition (plasmid-vs-endpoint is not a replicate pair).
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === CONFIGURE FOR YOUR SCREEN ===
# stage_map: sample -> 'plasmid' | 'day_0' | 'endpoint'   (drives the threshold band)
# condition_map: condition -> [replicate sample columns]  (drives which pairs are correlated)
stage_map = {'Plasmid': 'plasmid', 'Day0_r1': 'day_0', 'Day0_r2': 'day_0',
             'Endpoint_r1': 'endpoint', 'Endpoint_r2': 'endpoint'}
condition_map = {'day_0': ['Day0_r1', 'Day0_r2'], 'endpoint': ['Endpoint_r1', 'Endpoint_r2']}

STAGE_THRESHOLDS = {                      # same numbers as SKILL.md's stage_specific_thresholds()
    'plasmid':  {'pct_zero_max': 0.5, 'gini_pass': 0.10, 'gini_fail': 0.20},
    'day_0':    {'pct_zero_max': 1.0, 'gini_pass': 0.12, 'gini_fail': 0.25},
    'endpoint': {'pct_zero_max': 5.0, 'gini_pass': 0.30, 'gini_fail': 0.55},
}

# Load counts (MAGeCK format)
counts = pd.read_csv('screen.count.txt', sep='\t', index_col=0)


def validate_counts(counts, stage_map):
    '''Fail with a clear message instead of a traceback halfway through the report.'''
    if 'Gene' not in counts.columns:
        raise ValueError("count table has no 'Gene' column (MAGeCK format: sgRNA, Gene, samples...)")
    matrix = counts.drop('Gene', axis=1)
    non_numeric = [c for c in matrix.columns if not pd.api.types.is_numeric_dtype(matrix[c])]
    if non_numeric:
        raise ValueError(f'non-numeric sample columns: {non_numeric}')
    if (matrix.values < 0).any():
        raise ValueError('negative values in the count table; these are not read counts')
    if counts.index.duplicated().any():
        raise ValueError(f'{counts.index.duplicated().sum()} duplicated sgRNA identifiers')
    missing_stage = [c for c in matrix.columns if c not in stage_map]
    if missing_stage:
        raise ValueError(f'no stage assigned for samples {missing_stage}; add them to stage_map')
    dead = [c for c in matrix.columns if matrix[c].sum() == 0]
    if dead:
        print(f'WARNING: samples with zero total reads (sequencing failure), excluded: {dead}')
        matrix = matrix.drop(columns=dead)
    return matrix


genes = counts['Gene']
count_matrix = validate_counts(counts, stage_map)

print('=== CRISPR Screen QC Report ===\n')
print(f'sgRNAs: {len(count_matrix)}')
print(f'Genes: {genes.nunique()}')
print(f'Samples: {list(count_matrix.columns)}\n')

# Zero counts - the acceptable fraction depends on stage: a plasmid pool should have almost none,
# an endpoint sample legitimately loses guides to selection.
print('--- Library Representation (stage-specific) ---')
zero_pct = (count_matrix == 0).sum() / len(count_matrix) * 100
for sample, pct in zero_pct.items():
    limit = STAGE_THRESHOLDS[stage_map[sample]]['pct_zero_max']
    status = 'PASS' if pct <= limit else 'WARN' if pct <= limit * 2 else 'FAIL'
    print(f'{sample} [{stage_map[sample]}]: {pct:.2f}% zero counts (limit {limit}%) [{status}]')


# Gini index
def gini_index(x):
    x = np.sort(x[x > 0].astype(float))
    if x.size == 0:                      # a fully failed sample: report, do not raise
        return np.nan
    n = x.size
    cumx = np.cumsum(x)
    return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n


# Gini measures count inequality; the pass band differs by stage (SKILL.md, "Gini Coefficient")
print('\n--- Read Distribution (Gini, stage-specific) ---')
for sample in count_matrix.columns:
    g = gini_index(count_matrix[sample].values)
    band = STAGE_THRESHOLDS[stage_map[sample]]
    if np.isnan(g):
        status = 'NO DATA'
    else:
        status = 'PASS' if g < band['gini_pass'] else 'WARN' if g < band['gini_fail'] else 'FAIL'
    print(f'{sample} [{stage_map[sample]}]: {g:.3f} (pass <{band["gini_pass"]}, fail >{band["gini_fail"]}) [{status}]')

# Replicate correlation - within a condition only. Correlating a plasmid sample against an
# endpoint sample measures selection, not reproducibility, and drags the average down.
print('\n--- Replicate Correlation (within condition) ---')
log_counts = np.log10(count_matrix + 1)
corr = log_counts.corr()
for condition, samples in condition_map.items():
    samples = [s for s in samples if s in count_matrix.columns]
    for i, c1 in enumerate(samples):
        for c2 in samples[i + 1:]:
            r = corr.loc[c1, c2]
            rho = log_counts[c1].corr(log_counts[c2], method='spearman')
            status = 'PASS' if r > 0.8 else 'WARN' if r > 0.6 else 'FAIL'
            print(f'{condition}: {c1} vs {c2}: Pearson={r:.3f} Spearman={rho:.3f} [{status}]')
    if len(samples) < 2:
        print(f'{condition}: only {len(samples)} sample(s); no replicate pair to check')

# Plot
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for sample in count_matrix.columns:
    axes[0].hist(np.log10(count_matrix[sample] + 1), bins=50, alpha=0.5, label=sample)
axes[0].set_xlabel('Log10(counts + 1)')
axes[0].set_ylabel('sgRNAs')
axes[0].legend()

import seaborn as sns
sns.heatmap(corr, annot=True, cmap='RdYlBu_r', vmin=0.5, vmax=1, ax=axes[1])
axes[1].set_title('Sample Correlation (all pairs; only within-condition pairs are replicates)')

plt.tight_layout()
plt.savefig('screen_qc.png', dpi=150)
print('\nQC plots saved to screen_qc.png')
