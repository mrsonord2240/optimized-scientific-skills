---
name: bio-proteomics-proteomics-qc
description: Quality control for bottom-up proteomics across three levels -- instrument/raw-signal (mass accuracy, RT/iRT fit, FWHM, TIC vs injection time, % MS2 identified), identification/run (missed cleavages, charge states, PTM handling artifacts, contaminants), and experiment/quantitative (replicate correlation on log2, CV on the linear scale, completeness, MNAR-vs-MCAR missingness, PCA/batch, TMT channel balance, DIA q-values). Frames QC as a control chart against a per-instrument rolling baseline, not fixed cutoffs, and mandates inspecting raw boxplots, per-sample ID counts, total signal, and contaminant removal BEFORE normalizing -- because median normalization erases loading failures. Use when assessing proteomics data quality, diagnosing outlier samples, or deciding which samples to exclude before differential testing. The statistical test itself is differential-abundance; normalization mechanics are quantification; DIA q-value internals are dia-analysis.
tool_type: mixed
primary_tool: pandas
license: MIT
---

## Version Compatibility

Reference examples tested with: pandas 2.2+, numpy 1.26+, scipy 1.12+, matplotlib 3.8+, scikit-learn 1.4+, limma 3.58+ (re-checked on pandas 3.0.5, numpy 2.5.3, scikit-learn 1.9.1, PTXQC 1.1.5 with rmarkdown 2.31 and Pandoc 3.11)

Seed anything stochastic and re-run it once before you report it: `sklearn.decomposition.PCA` has no default `random_state` and its `'auto'` solver is randomized on wide matrices, so an unseeded QC PCA prints different p-values on identical data.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Proteomics Quality Control -- A Three-Level Funnel Where the Matrix Sees the Failure Last

**"Check the quality of my proteomics data"** -> Read instrument, identification, and quantitative metrics as a descending funnel of silent failures, and inspect raw signal BEFORE normalizing -- because by the time a fault reaches the deliverable matrix, normalization has usually erased the evidence.
- Python: `pandas` for matrix QC; `matplotlib`/`seaborn` for raw boxplots, correlation heatmaps, PCA
- R: PTXQC `createReport()` for MaxQuant search-table QC; `limma::plotMDS()`/`plotDensities()`; MSstatsTMT `dataProcessPlotsTMT()` for TMT channel balance

Scope: This skill OWNS QC diagnosis across all three levels -- which metric localizes which fault, what threshold means trouble, and the mandatory inspect-before-normalize ordering. Normalization mechanics route to quantification; the differential test routes to differential-abundance; DIA q-value computation routes to dia-analysis. OUT OF SCOPE: running the statistical test, the normalization algorithms themselves, and DIA q-value/FDR internals.

## The Single Most Important Modern Insight -- QC Is a Three-Level Funnel and Normalization Hides the Evidence

1. **QC is a three-level funnel of silent failures, and the deliverable matrix is the LAST place a problem becomes visible.** Faults originate at the instrument (spray, calibration, column) or in identification (digestion, contamination, PTM artifacts), but a protein matrix only shows the downstream symptom -- a low correlation or an outlier sample. A matrix-only QC pass is one-third of the job and blind to where faults actually start. Localize by descending: read every metric together with its co-readouts, never alone.

2. **Almost no metric has a universal pass/fail cutoff; the defensible practice is a per-instrument, per-method control chart.** Deviation from a lab's own rolling baseline (Levey-Jennings, +/-2 SD warn, +/-3 SD action) detects faults that a constant threshold misses or false-flags (Neely and Palmblad 2024). The numbers below seed a control chart, they are not standards-body limits.

3. **Median normalization HIDES loading problems that must be SEEN first.** Median/quantile normalization works by forcing a chosen summary statistic of every sample equal. A sample that genuinely loaded 3x low sits visibly shifted down in a RAW boxplot -- an obvious, diagnosable defect. The instant the matrix is median-normalized, the algorithm shifts that sample up by a constant to match everyone's median; boxplots line up perfectly; the evidence is mathematically erased. Worse, the low-loaded sample's noisy low-abundance signal gets stretched up to mid-range and injected into the differential test while QC plots look pristine. MANDATE: inspect raw/un-normalized boxplots plus per-sample ID counts, total signal, and missing fraction BEFORE normalizing; remove loading/injection failures and contaminants; THEN normalize and re-plot on the survivors. The identical principle governs TMT channel-loading balance.

## The Three Levels in One Table

| Level | Question | Inputs | Faults localized |
|-------|----------|--------|------------------|
| 1. Instrument / raw-signal | Is the LC-MS hardware performing? | Vendor `.raw`/`.d`; RawTools/RawBeans/rawrr/rawDiag/QuaMeter; Panorama AutoQC | Column, spray/emitter, mass analyzer/calibration |
| 2. Identification / run | Did this run identify peptides correctly? | Search tables (MaxQuant `txt/`, FragPipe `*.tsv`, DIA-NN report); PTXQC | Digestion, sample-handling PTM artifacts, contamination, FDR efficiency |
| 3. Experiment / quantitative | Are the numbers reproducible and comparable? | Protein/peptide intensity matrix; MSstatsTMT, pandas/limma | Loading/pipetting, batch, outliers, missingness, sample swaps |

The most integrative metric (% MS2 identified / ID count) is the first alarm but the LEAST specific -- it moves whenever anything upstream degrades. The same protein-count drop means spray (erratic TIC + maxed injection time), column (lost RT + broad peaks + rising backpressure), or sample (high contaminant fraction) depending on what co-moves.

## Tool Taxonomy

| Tool / method | Level | Citation | Mechanism / role | When |
|---------------|-------|----------|------------------|------|
| PTXQC (R, CRAN) | 2+3 | Bielow 2016 | `createReport()` over MaxQuant `txt/` or mzTab; per-metric scores in [0,1], QC heatmap PDF (+ mzQC + YAML; HTML only when Pandoc is reachable) | MaxQuant output, fast multi-metric report |
| RawTools / RawBeans | 1 | Kovalchik 2019; Morgenstern 2021 | Parse Thermo `.raw` for IT, TIC, FWHM, scan timing | Diagnose instrument faults from raw files |
| rawrr / rawDiag | 1 | Kockmann 2021; Trachsel 2018 | R access to Orbitrap scan metadata | Custom Level-1 plots / method optimization |
| QuaMeter | 1+2 | Ma 2012 | Vendor-independent ID-free and ID-based metrics | Cross-vendor Level-1 QC |
| Skyline + Panorama AutoQC | 1 (longitudinal) | Bereman 2016 | Levey-Jennings + CUSUM/Moving-Range, SD-band flagging | System-suitability trending over time |
| MSstatsTMT | 3 (TMT) | Huang 2020 | `proteinSummarization()`, `dataProcessPlotsTMT()`; filters isolation interference on import | TMT channel-balance and QC plots |
| pandas / limma matrix QC | 3 | (this skill) | Correlation, CV, completeness, PCA on the matrix | Experiment-level QC (the code below) |
| differential-abundance | (route OUT) | -- | The moderated test itself | Hit calling after QC passes |
| quantification | (route OUT) | -- | Normalization and imputation mechanics | The how of normalizing |
| dia-analysis | (route OUT) | -- | DIA q-value/FDR internals | DIA-NN/Spectronaut report computation |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| MaxQuant `txt/` folder, want fast multi-metric report | PTXQC `createReport(txt_folder=...)` | Scores Level-2/3 metrics vs a representative file; writes the PDF + mzQC + heatmap + YAML with no extra dependency, and an HTML report ONLY when Pandoc is reachable |
| One run per condition (no replicates) | Report which checks could NOT run: replicate r, CV and the within-group loading rule all need >=2 samples per group | An empty correlation table and a NaN CV mean "not measured", never "clean" |
| Protein-count drop, cause unknown | Descend to Level 1: read TIC + injection time + RT/FWHM together | Co-readouts localize spray vs column vs sample |
| Replicate correlation low for one sample | Check if it correlates better with a DIFFERENT group | Distinguishes sample swap from prep failure |
| Boxplots flat but a sample feels wrong | Re-plot the RAW (un-normalized) matrix | Normalization erased the loading evidence |
| Deciding how to impute | Diagnose MNAR (left tail) vs MCAR (all-abundance) from the histogram FIRST | Wrong imputer corrupts present/absent calls |
| TMT data, channel looks off | Within-plex channel totals on raw reporter intensities (code below); MSstatsTMT `dataProcessPlotsTMT` only after `proteinSummarization(..., global_norm = FALSE, reference_norm = FALSE)` | The MSstatsTMT default (`global_norm = TRUE`) equalizes channel medians, so its QC plot hides imbalance; it also needs PSM-level input |
| DIA matrix, how many proteins are real | Filter Global.Q.Value and Global.PG.Q.Value, route q internals to dia-analysis | Precursor q != protein q; both needed |
| Long sample queue, drift suspected | Interspersed QC every 4th-5th injection + Levey-Jennings | Turns one check into a time series |

Default when uncertain: plot the RAW per-sample boxplots, ID counts, total signal, and missing fraction first; remove loading/injection failures and contaminants; only then normalize, re-plot, and proceed to correlation/CV/PCA on the survivors.

## Inspect Raw Signal and Remove Contaminants Before Normalizing

**Goal:** Catch loading/injection failures and strip contaminant/decoy rows while they are still visible -- before normalization erases them.

**Approach:** Load the un-normalized matrix, plot per-sample boxplots plus ID counts and total signal, filter MaxQuant `Potential contaminant`/`Reverse`/`Only identified by site` rows, THEN log-transform and normalize on the survivors.

Which column is un-normalized: MaxQuant `Intensity <sample>` (`LFQ intensity <sample>` is already MaxLFQ-normalized and hides a 3x-low load); DIA-NN `Precursor.Quantity` (not `Precursor.Normalised` or `PG.MaxLFQ`); TMT raw reporter intensities before any channel/global normalization. Search engines write missing values as 0 (MaxQuant intensities, DIA-NN `PG.MaxLFQ`), so convert 0 to NaN before counting IDs or missingness.

```python
import pandas as pd
import numpy as np

contaminant_flags = ['Potential contaminant', 'Reverse', 'Only identified by site']

def strip_contaminant_rows(protein_groups):
    keep = pd.Series(True, index=protein_groups.index)
    for col in contaminant_flags:
        match = next((c for c in protein_groups.columns if c.lower() == col.lower()), None)  # MaxQuant casing varies by version -- match case-insensitively
        if match is not None:
            keep &= protein_groups[match].fillna('') != '+'  # MaxQuant marks flagged rows with a literal '+'
    return protein_groups[keep]

def raw_sample_qc(raw_intensities, sample_groups):
    raw = raw_intensities.replace(0, np.nan)  # MaxQuant/DIA-NN write missing as 0; zeros are NOT quantified
    qc = pd.DataFrame({
        'n_quantified': raw.notna().sum(),
        'total_signal': raw.sum(),
        'median_intensity': raw.median(),
        'missing_pct': 100 * raw.isna().sum() / len(raw)})
    group = sample_groups.reindex(qc.index)
    sizes = group.value_counts()
    singletons = sorted(sizes[sizes < 2].index.astype(str))
    if singletons:
        # A group of one IS its own median: fold and ids would be exactly 1.000 and the flag
        # could never fire. Say so and fall back to the all-sample baseline for those samples.
        print(f'WARNING: single-sample group(s) {singletons}: the within-group loading rule cannot '
              'fire there, so those samples are compared to the ALL-sample median instead. A loading '
              'difference that tracks condition is NOT detectable in a design with no replicates.')
        baseline = group.where(group.map(sizes) >= 2, 'ALL')
    else:
        baseline = group
    qc['baseline'] = baseline
    qc['fold_total_vs_group'] = qc['total_signal'] / qc.groupby(baseline)['total_signal'].transform('median')
    qc['ids_vs_group'] = qc['n_quantified'] / qc.groupby(baseline)['n_quantified'].transform('median')
    qc['flag'] = (qc['fold_total_vs_group'] <= 0.5) | (qc['ids_vs_group'] < 0.8)  # >=2x low total, or >20% fewer IDs
    return qc

def contaminant_fraction(protein_groups, intensity_cols, flag_col='Potential contaminant'):
    flagged = protein_groups[flag_col].fillna('') == '+'
    raw = protein_groups[intensity_cols].replace(0, np.nan)
    return 100 * raw[flagged].sum() / raw.sum()  # percent of summed raw intensity, per sample
```

Read the boxplots before normalizing, but apply the loading rule to TOTAL raw signal and ID count, not the boxplot median: a sample with total signal >=2x below its group median, or an ID count more than 15-20% below it, is a loading/injection failure to exclude, not to rescale. Left-censoring removes a low-loaded sample's weakest values, so its median looks less shifted than it is (synthetic test: 0.41x total but 0.62x median). Judge the contaminant fraction against the lab's own baseline for that sample type and check whether it differs between groups; there is no universal cutoff (PTXQC's 1% threshold belongs to its user-defined special-contaminant plot, not the general contaminant score). Keratin and trypsin autolysis dominate LOW-INPUT samples (single-cell, IPs, gel bands) because they are a roughly fixed absolute amount whose fractional share explodes as load shrinks.

With ONE run per condition the within-group form of this rule is inert -- every sample is its own group median, so `fold_total_vs_group` is exactly 1.0 and nothing can flag. The function above falls back to the all-sample median and prints a warning; report that fallback, and report that a loading difference which tracks condition is not separable from biology in that design.

## Replicate Correlation on log2

**Goal:** Quantify reproducibility without letting a few abundant proteins fake agreement.

**Approach:** Correlate on log2 intensities (variance-stabilized, high-abundance tail compressed), report within-group pairs, and flag a sample correlating better with another group as a possible swap.

```python
from itertools import combinations

def replicate_correlation(log2_intensities, sample_groups):
    # DESIGN REQUIREMENT: >=2 samples in at least one group. Without it this silently returned
    # an empty table, which reads like "no reproducibility problem" when it means "not measured".
    sizes = sample_groups.value_counts()
    if (sizes >= 2).sum() == 0:
        raise ValueError(f'replicate_correlation needs >=2 samples in a group; sizes are '
                         f'{sizes.to_dict()}. With one run per condition there is no replicate '
                         'reproducibility to measure -- do not report "no outliers found".')
    unchecked = sorted(sizes[sizes < 2].index.astype(str))
    if unchecked:
        print(f'WARNING: no within-group pairs for {unchecked}; those samples are UNCHECKED here.')
    corr = log2_intensities.corr(method='pearson')  # log2 first: Pearson on raw is a high-abundance artifact
    rows = []
    for group in sample_groups.unique():
        members = sample_groups[sample_groups == group].index
        for s1, s2 in combinations(members, 2):
            rows.append({'group': group, 's1': s1, 's2': s2, 'r': corr.loc[s1, s2]})
    return pd.DataFrame(rows)
```

Technical replicates r > 0.98 (instrument noise only); biological r ~ 0.90-0.98 (genuine variance, lower is expected and correct); soft floor r > 0.8 to retain a biological replicate. A Spearman check is a robustness aid only -- ranks discard the magnitude that quant QC cares about.

## Coefficient of Variation on the Linear Scale

**Goal:** Summarize per-condition precision with a number that means what it says.

**Approach:** Compute CV = SD/mean on LINEAR (non-log) intensities; if only logged values exist use the geometric-CV formula. Report the median CV per condition (the per-protein distribution is right-skewed).

```python
def median_cv_linear(linear_intensities, sample_groups):
    # DESIGN REQUIREMENT: >=2 samples per group. A one-sample group yields NaN, and a table of
    # NaNs is not "excellent precision" -- fail loudly instead of returning it.
    sizes = sample_groups.value_counts()
    if (sizes >= 2).sum() == 0:
        raise ValueError(f'median_cv_linear needs >=2 samples in a group; sizes are {sizes.to_dict()}. '
                         'CV is undefined with no replicates -- report "not measurable", not NaN.')
    rows = []
    for group in sample_groups.unique():
        members = sample_groups[sample_groups == group].index
        if len(members) < 2:
            print(f'WARNING: group {group!r} has {len(members)} sample(s); its CV is NaN (undefined), not low.')
        block = linear_intensities[members]
        per_protein_cv = block.std(axis=1) / block.mean(axis=1)  # base CV formula REQUIRES linear scale
        rows.append({'group': group, 'median_cv_pct': 100 * per_protein_cv.median()})
    return pd.DataFrame(rows)

def geometric_cv_from_log(log_intensities):
    sigma = log_intensities.std(axis=1) * np.log(2)  # convert log2 SD to natural-log SD
    return 100 * np.sqrt(np.expm1(sigma ** 2))  # gCV = sqrt(exp(sigma^2) - 1)
```

Applying the base formula to log-transformed data compresses CV by roughly ln2 x mean log2 intensity -- about 15-20x at typical MaxQuant intensities, scale-dependent (most proteins appear to have CV < 1%) -- meaningless (Brenes 2024). State normalization state, transform, and software params or the CV is uninterpretable: DIA-NN "High precision" mode silently median-normalizes, halving median CV vs "High accuracy". Technical median CV < ~10-20%, biological ~20-40%; a LOWER CV is not automatically better (loose FDR or faulty MS1 extraction produce artificially low CVs).

## Missingness Mechanism and Completeness

**Goal:** Decide how to impute by first deciding why values are missing.

**Approach:** Diagnose the missingness profile -- left-tail concentration means MNAR (left-censored, abundance-dependent), all-abundance scatter means MCAR -- and filter on completeness before imputing only the shallow remainder.

```python
def missingness_profile(log2_intensities, n_bins=10):
    present_per_protein = log2_intensities.notna().mean(axis=1)
    mean_abundance = log2_intensities.mean(axis=1)
    abundance_bin = pd.qcut(mean_abundance, n_bins, duplicates='drop')
    # present fraction per mean-abundance bin: rising-with-abundance = MNAR, flat = MCAR
    return present_per_protein.groupby(abundance_bin, observed=True).mean()

def completeness_filter(log2_intensities, sample_groups, min_valid_frac=0.7):
    keep = pd.Series(False, index=log2_intensities.index)
    for group in sample_groups.unique():
        block = log2_intensities[sample_groups[sample_groups == group].index]
        keep |= block.notna().mean(axis=1) >= min_valid_frac  # valid in >=70% of >=1 condition
    return log2_intensities[keep]
```

kNN-imputing a genuinely-absent (MNAR) value invents mid-range abundance and KILLS a real present/absent difference; a left-shifted draw (Perseus down-shifted normal, downshift=1.8 SD below the observed mean, width=0.3 of observed SD) on an MCAR gap FABRICATES a false low and inflates a difference. Match the imputer to the mechanism. The imputation mechanics themselves are quantification.

## PCA and Batch Detection

**Goal:** See whether the dominant variance is biology or batch, and flag outlier samples.

**Approach:** On the normalized survivors, run PCA, color by condition and by batch, and test whether top PCs associate with batch.

```python
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import f_oneway

def pca_batch_check(normalized_log2, sample_info, batch_col='batch'):
    # sample_info must be indexed by sample name, e.g. pd.read_csv(...).set_index('sample')
    if set(sample_info.index) != set(normalized_log2.columns):
        raise ValueError('sample_info index must equal the matrix column names (set_index on the sample column)')
    # complete cases only: a row-median fill pulls high-missing (failed) samples to the centre of the PCA
    complete = normalized_log2.dropna(how='any')
    n_samples = complete.shape[1]
    if n_samples < 3 or len(complete) < n_samples:
        raise ValueError(f'too few samples ({n_samples}) or complete proteins ({len(complete)}) for PCA')
    n_pc = min(5, n_samples - 1)
    scaled = StandardScaler().fit_transform(complete.T)
    # REPRODUCIBILITY: sklearn's default svd_solver='auto' switches to RANDOMIZED SVD on a wide
    # matrix (samples x proteins), so an unseeded PCA returns different minor components and a
    # different printed p-value on every identical run. 'full' is exact and cheap at QC sizes;
    # random_state pins the randomized path for anyone who changes the solver back.
    pcs = PCA(n_components=n_pc, svd_solver='full', random_state=0).fit(scaled)
    coords = pd.DataFrame(pcs.transform(scaled), columns=[f'PC{i+1}' for i in range(n_pc)],
                          index=complete.columns).join(sample_info)
    print(f'PCA on {len(complete)} complete proteins of {len(normalized_log2)}')
    for pc in coords.columns[:min(3, n_pc)]:
        groups = [coords[coords[batch_col] == b][pc] for b in coords[batch_col].unique()]
        # a level with one sample makes f_oneway raise 'At least two samples are required; got 1'
        if len(groups) < 2 or min(len(g) for g in groups) < 2:
            print(f'{pc} ~ {batch_col}: NOT TESTABLE, level sizes {[len(g) for g in groups]} '
                  f'(need >=2 levels with >=2 samples each) -- this is not evidence of no batch effect')
            continue
        _, p = f_oneway(*groups)
        print(f'{pc} ~ {batch_col}: p={p:.4f}')
    return coords, pcs.explained_variance_ratio_
```

A sample isolated from its group is a removal/re-run candidate, but a high-missing sample is judged by `raw_sample_qc`, not by PCA. If batch is PC1, keep batch in the design matrix for the differential test (preferred when batch and condition are balanced), and use `limma::removeBatchEffect` (or ComBat) only on the matrix used for PCA/plots to re-inspect biology; do not test on a batch-corrected matrix and also model batch. If batch is FULLY confounded with condition (every batch level holds exactly one condition) nothing can be corrected: batch and condition are the same variable, and removing one removes the other -- on a fully confounded synthetic set the mean |log2FC| of 104 truly-changed proteins went from 1.55 to 0.00 after batch removal. Report the design as non-identifiable and stop; do not correct, and do not test. Stop and ask before excluding samples, when n < 5 per group makes PCA unstable, or when no un-normalized column is available for the loading check. Visualization of the projection routes to data-visualization/dimensionality-reduction-plots.

## TMT Channel Balance Within Each Plex

**Goal:** Catch an under- or over-loaded TMT channel before any channel normalization.

**Approach:** Per plex, compare each channel's summed raw reporter intensity to the plex median; plexes differ in overall signal, so never compare channels across plexes directly (use the reference-channel ratio for that).

```python
def tmt_channel_balance(plex_matrices):
    # plex_matrices: {plex_name: DataFrame proteins x channels of RAW reporter intensities}
    rows = []
    for plex, m in plex_matrices.items():
        total = m.replace(0, np.nan).sum()
        for channel, fold in (total / total.median()).items():
            rows.append({'plex': plex, 'channel': channel, 'fold_vs_plex_median': fold})
    balance = pd.DataFrame(rows)
    balance['investigate'] = np.abs(np.log2(balance['fold_vs_plex_median'])) > 1  # > 2x; flag > 3-4x
    return balance
```

## Per-Method Failure Modes

### Median normalization hides loading failures
**Trigger:** Normalizing the matrix before inspecting raw per-sample signal, or "inspecting" an already-normalized column (MaxQuant `LFQ intensity`, DIA-NN `PG.MaxLFQ`/`Precursor.Normalised`). **Mechanism:** median-centering shifts each sample by a constant to equalize the very statistic that was the symptom of a low load. **Symptom:** flat, clean boxplots that hide a 3x-low sample now stretched into mid-range. **Fix:** plot RAW boxplots + ID counts + total signal first (MaxQuant `Intensity`, DIA-NN `Precursor.Quantity`); exclude failures; then normalize.

### Normalizing with contaminants still in the matrix
**Trigger:** Contaminant/decoy rows left in before log + normalize. **Mechanism:** keratin/trypsin/albumin inflate the denominator and shift the median; when their load differs across groups the differential gets normalized into the real proteins. **Symptom:** spurious fold changes; a contaminant fraction that varies by group. **Fix:** filter `Potential contaminant` + `Reverse` + `Only identified by site` BEFORE log + normalize.

### Wrong imputer for the missingness mechanism
**Trigger:** kNN on MNAR, or left-shift on MCAR. **Mechanism:** kNN borrows mid-range neighbors for a value that is low because it is absent; left-shift draws a deep low for a value missing at random. **Symptom:** killed present/absent calls (kNN-on-MNAR) or inflated false lows (left-shift-on-MCAR). **Fix:** diagnose left-tail vs all-abundance from the histogram first; mechanics route to quantification.

### CV computed on log-transformed data
**Trigger:** Base CV formula applied after log2. **Mechanism:** SD/mean is defined for linear intensity; logging compresses it by about ln2 x mean log2 intensity (15-20x at typical MaxQuant intensities). **Symptom:** most proteins appear to have CV < 1%. **Fix:** compute on linear intensity, or use the geometric-CV formula on logged values; always state transform + normalization + software.

### Pearson r on raw intensity
**Trigger:** Correlating un-logged intensities. **Mechanism:** a few high-abundance proteins dominate the covariance. **Symptom:** r = 0.99 while the bulk disagrees. **Fix:** log2 before correlating; Spearman as a robustness check only.

### Bimodal mass-error histogram read as calibration
**Trigger:** A two-peaked ppm-error distribution. **Mechanism:** almost always monoisotopic mis-assignment or co-isolation (a search/sample problem), NOT calibration drift; a generous tolerance still IDs the mis-assigned precursors so ID rate looks fine. **Symptom:** bimodal histogram, normal ID rate. **Fix:** correct isotope-error tolerance/deisotoping, not recalibration.

### Charge-state distribution off the platform baseline
**Trigger:** The fully-tryptic 2+ fraction drifts from the rolling baseline. **Mechanism:** a tryptic peptide carries two basic sites (C-terminal K/R + N-terminus) so 2+ dominates; excess 3+/4+ comes from internal basic residues left by missed cleavages, excess 1+ from poor ionization, short peptides, or contaminants. **Symptom:** raised high-charge fraction (a digestion/chemistry signal) or raised 1+ fraction (an ionization/spray signal). **Fix:** read the charge distribution together with the missed-cleavage rate (high charge co-moving with missed cleavages = under-digestion) to separate a chemistry problem from a spray problem a raw protein-count drop cannot resolve alone.

### TMT ratio compression from co-isolation
**Trigger:** Isobaric quant with a wide isolation window. **Mechanism:** near-isobaric co-eluting precursors are co-isolated and add their own reporters across all channels, a uniform pedestal. **Symptom:** every fold-change compressed toward 1 (a real 10:1 reads ~5:1). **Fix:** filter isolation interference < 50%; prefer SPS-MS3 (McAlister 2014) / FAIMS / narrow windows; run a TKO or empty-channel control to measure the floor.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Mass accuracy (internal/lock) median \|err\| < 1-3 ppm, single-mode centered 0 | CONVENTION | width approaching MS1 tolerance loses real IDs |
| iRT RT-fit R^2 > 0.99 (warn below) | CONVENTION (mechanism firm) | residual is just LC noise on a stable gradient |
| FWHM alarm on > 20-30% rise vs baseline; >= 8-10 points across FWHM (floor ~5) | Kocher 2011 | peak capacity tracks peptide IDs; points needed for accurate AUC |
| % MS2 identified: < 20% bad, 20-35% ok, >= 35% great | PTXQC `createYaml.R` | generic pass marks, not a biological ceiling |
| Missed cleavages: >= 75-85% at 0 MC; flag > 25-30% with >= 1 MC | OPERATIONAL (PTXQC-style) | porcine trypsin ~78% efficient even ideally |
| Replicate Pearson r (log2): technical > 0.98, biological 0.90-0.98, floor 0.8 | CONVENTION (mechanism firm) | log2 variance-stabilizes; biological variance is real |
| Median CV (linear): technical < 10-20%, biological 20-40% | CONVENTION | DIA < DDA (no stochastic sampling); lower not always better |
| Completeness: valid in >= 50-70% of replicates in >= 1 condition | CONVENTION | filter before imputing |
| Perseus MNAR imputation: downshift = 1.8 SD, width = 0.3 | Tyanova 2016 | deep left tail simulates below-LOD, narrowed so not mistaken for real |
| TMT channel deviation: investigate > ~2x, flag > ~3-4x | CONVENTION | pipetting/labeling vs biology |
| Isolation interference < 50% PSM filter (< 30% stricter) | CONVENTION (PD practice) | above ~50% contaminant dominates, ratios uninterpretable |
| DIA precursor + protein q both <= 0.01, GLOBAL q, PICKED protein estimator | STANDARD | precursor != protein FDR; route internals to dia-analysis |
| Levey-Jennings: +/-2 SD warn, +/-3 SD action; QC every 4th-5th injection | Bereman 2016 | ~95/99.7% of points under stable normal; catch drift early |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `KeyError: 'Mass Error [ppm]'` | MaxQuant column casing varies by version | match case-insensitively (`Mass error` vs `Mass Error`); never hard-code |
| Contaminant rows have `True`/`False`, filter keeps all | MaxQuant flags with a literal `'+'`, not a boolean | filter `col != '+'` |
| CV unexpectedly tiny (< 1%) | base CV formula applied to log2 data | compute on linear intensity or use geometric CV |
| `r = 0.99` but samples clearly differ | Pearson on raw (un-logged) intensity | log2 transform before correlating |
| PCA dominated by injection day | batch effect, not biology | batch in the design for the test; `removeBatchEffect`/ComBat only for plots, then re-inspect |
| Every sample shows 0% missing and identical ID counts | MaxQuant/DIA-NN zeros counted as values | `.replace(0, np.nan)` before counting |
| `TypeError: At least two samples are required; got 1` in `pca_batch_check` | Either `sample_info` has a RangeIndex so the join produced NaN batches, OR a level of `batch_col` genuinely has one sample (e.g. one run per condition) | `sample_info.set_index('sample')` first. If the index is already correct the design is the cause: the function now prints `NOT TESTABLE` with the level sizes instead of raising -- report that, do not read it as "no batch effect" |
| `pca_batch_check` prints a different p-value on every identical run | sklearn's default `svd_solver='auto'` picks randomized SVD on a wide matrix, and `PCA()` has no `random_state` | `PCA(n_components=n_pc, svd_solver='full', random_state=0)`; PC1/PC2 were already stable, PC3+ and every printed p-value were not (six unseeded fits on a 20x600 matrix gave six different PC3 p-values) |
| `replicate_correlation` returns 0 rows / `median_cv_linear` returns NaN | Fewer than 2 samples in every group | Not a clean result -- both now raise with the group sizes. Say which checks could not run |
| PTXQC prints "The 'Pandoc' converter is not installed on your system ... Pandoc is required for HTML reports" | `createReport()` renders the HTML leg through rmarkdown, which needs Pandoc on `PATH` or at `RSTUDIO_PANDOC` | NOT fatal: the PDF, mzQC, heatmap and YAML are still written and the call completes. For the HTML too, install Pandoc or `Sys.setenv(RSTUDIO_PANDOC='<pandoc dir>')` and confirm `rmarkdown::pandoc_available()` is TRUE before the call |
| MSstatsTMT QC plot shows identical channel medians | `proteinSummarization` default `global_norm = TRUE` | re-run with `global_norm = FALSE, reference_norm = FALSE` for the balance view |
| PTXQC "not found" via `BiocManager` | PTXQC is on CRAN, not Bioconductor | `install.packages('PTXQC')` |
| `createReport()` errors on a dataframe arg | it takes a txt-folder path / mzTab / YAML, not dataframes | pass `txt_folder=` (the MaxQuant `txt/` directory) |

## References

- Bielow C, Mastrobuoni G, Kempa S. Proteomics Quality Control: Quality Control Software for MaxQuant Results. *J Proteome Res* 2016;15(3):777-787.
- Kovalchik KA, Colborne S, Spencer SEP, et al. RawTools: Rapid and Dynamic Interrogation of Orbitrap Data Files for Mass Spectrometer System Management. *J Proteome Res* 2019;18(2):700-708.
- Morgenstern D, Barzilay R, Levin Y. RawBeans: A Simple, Vendor-Independent, Raw-Data Quality-Control Tool. *J Proteome Res* 2021;20(4):2098-2104.
- Kockmann T, Panse C. The rawrr R Package: Direct Access to Orbitrap Data and Beyond. *J Proteome Res* 2021;20(4):2028-2034.
- Trachsel C, Panse C, Kockmann T, et al. rawDiag: An R Package Supporting Rational LC-MS Method Optimization for Bottom-up Proteomics. *J Proteome Res* 2018;17(8):2908-2914.
- Ma ZQ, Polzin KO, Dasari S, et al. QuaMeter: Multivendor Performance Metrics for LC-MS/MS Proteomics Instrumentation. *Anal Chem* 2012;84(14):5845-5850.
- Bereman MS, Beri J, Sharma V, et al. An Automated Pipeline to Monitor System Performance in Liquid Chromatography-Tandem Mass Spectrometry Proteomic Experiments. *J Proteome Res* 2016;15(12):4763-4769.
- Kocher T, Swart R, Mechtler K. Ultra-High-Pressure RPLC Hyphenated to an LTQ-Orbitrap Velos Reveals a Linear Relation between Peak Capacity and Number of Identified Peptides. *Anal Chem* 2011;83(7):2699-2704.
- Tyanova S, Temu T, Sinitcyn P, et al. The Perseus computational platform for comprehensive analysis of (prote)omics data. *Nat Methods* 2016;13(9):731-740.
- Brenes AJ. Calculating and Reporting Coefficients of Variation for DIA-Based Proteomics. *J Proteome Res* 2024;23(12):5274-5278.
- McAlister GC, Nusinow DP, Jedrychowski MP, et al. MultiNotch MS3 Enables Accurate, Sensitive, and Multiplexed Detection of Differential Expression across Cancer Cell Line Proteomes. *Anal Chem* 2014;86(14):7150-7158.
- Huang T, Choi M, Tzouros M, et al. MSstatsTMT: Statistical Detection of Differentially Abundant Proteins in Experiments with Isobaric Labeling and Multiple Mixtures. *Mol Cell Proteomics* 2020;19(10):1706-1723.
- Neely BA, Palmblad M, et al. Quality Control in the Mass Spectrometry Proteomics Core: A Practical Primer. *J Biomol Tech* 2024;35(3).

## Related Skills

- data-import - Load search-engine output and intensity matrices before QC
- quantification - Normalization and imputation mechanics that QC mandates running AFTER inspection
- differential-abundance - The moderated statistical test QC gates
- dia-analysis - DIA q-value/FDR internals behind the protein-count QC
- data-visualization/dimensionality-reduction-plots - PCA/MDS projection plotting
- workflows/proteomics-pipeline - End-to-end pipeline placing QC before differential testing
