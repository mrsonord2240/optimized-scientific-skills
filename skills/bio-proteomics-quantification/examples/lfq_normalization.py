'''Label-free normalization and TMT cross-plex IRS bridge, self-contained.'''
# Reference: numpy 2.5.3, pandas 3.0.5 | Verify API if version differs
import numpy as np
import pandas as pd

rng = np.random.default_rng(0)


def median_center(log_int):
    '''Subtract each sample median log2 intensity (corrects LOCATION only, cannot manufacture variance).'''
    sample_medians = log_int.median(axis=0)
    return log_int - sample_medians + sample_medians.median()


def sample_loading_normalize(plex):
    '''Per-channel scalar correcting total load WITHIN one plex; target = mean column sum.'''
    target = plex.sum(axis=0).mean()
    return plex * (target / plex.sum(axis=0))


def irs_scale(plexes, ref_cols):
    '''Internal Reference Scaling bridge (Plubell 2017): pin each plex reference channel to a common per-protein value.'''
    refs = pd.concat([p[ref] for p, ref in zip(plexes, ref_cols)], axis=1)
    refs = refs.where(refs > 0)    # a 0 or missing reference cannot anchor the bridge
    unbridged = refs.index[refs.isna().any(axis=1)]
    if len(unbridged):
        print(f'IRS: {len(unbridged)} proteins lack a reference in >=1 plex; set to NaN: {list(unbridged[:10])}')
    geomean = np.exp(np.log(refs).mean(axis=1, skipna=False))    # per-protein geometric mean across plexes
    out = []
    for i, p in enumerate(plexes):
        factor = geomean / refs.iloc[:, i]    # per-protein per-plex scaling factor
        out.append(p.mul(factor, axis=0))
    return out


def plex_offset(a, b, ref_a, ref_b):
    '''Median per-protein log2 offset between plexes on NON-reference channels (IRS equalizes references by construction).'''
    la = np.log2(a.drop(columns=ref_a)).mean(axis=1)
    lb = np.log2(b.drop(columns=ref_b)).mean(axis=1)
    return float(np.nanmedian(lb - la))


proteins = [f'P{i:03d}' for i in range(200)]
samples = [f'S{i}' for i in range(6)]

# Label-free: raw intensities with per-sample loading offsets and ~15% MaxQuant zeros (missing)
true_abundance = rng.normal(20, 2, (200, 6))
loading_offset = np.array([0.0, 0.4, -0.3, 0.6, -0.5, 0.2])    # per-sample load differences to be normalized out
raw = 2 ** (true_abundance + loading_offset)
raw[rng.random((200, 6)) < 0.15] = 0    # MaxQuant writes 0 for "not quantified"
lfq = pd.DataFrame(raw, index=proteins, columns=samples)

log_int = np.log2(lfq.replace(0, np.nan))    # 0 -> NaN before transform; log2(0) = -inf otherwise
print(f'Missing values: {100 * log_int.isna().sum().sum() / log_int.size:.1f}%')
print(f'Sample medians before centering: {np.round(log_int.median().values, 2)}')

normalized = median_center(log_int)
print(f'Sample medians after centering:  {np.round(normalized.median().values, 2)}')

# Filter proteins quantified in too few samples; min-peptides logic mirrors min-ratio-count=2 robustness
valid_per_protein = normalized.notna().sum(axis=1)
min_valid = len(samples) // 2    # require presence in >=50% of samples
filtered = normalized[valid_per_protein >= min_valid]
print(f'Proteins after >=50% valid filter: {len(filtered)}')

# TMT cross-plex: two 4-channel plexes over the same proteins, last channel a pooled reference.
# Plex B sees each protein at a random per-protein offset around 2x (elution sampling, not biology).
protein_level = rng.normal(18, 1.5, (200, 1))
plex_offset_true = rng.normal(1.0, 0.5, (200, 1))
plex_a = pd.DataFrame(2 ** (protein_level + rng.normal(0, 0.15, (200, 4))), index=proteins, columns=['A1', 'A2', 'A3', 'A_ref'])
plex_b = pd.DataFrame(2 ** (protein_level + plex_offset_true + rng.normal(0, 0.15, (200, 4))), index=proteins, columns=['B1', 'B2', 'B3', 'B_ref'])
plex_b.iloc[0, 3] = 0    # one protein without a plex-B reference reporter: reported and left unbridged
sl = [sample_loading_normalize(plex_a), sample_loading_normalize(plex_b)]
bridged = irs_scale(sl, ['A_ref', 'B_ref'])
print(f'Non-reference plex offset (median log2 B-A): SL only {plex_offset(sl[0], sl[1], "A_ref", "B_ref"):+.3f}, '
      f'after IRS {plex_offset(bridged[0], bridged[1], "A_ref", "B_ref"):+.3f}')
