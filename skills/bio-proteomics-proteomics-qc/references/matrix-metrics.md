# Matrix-level QC metrics: replicate correlation, CV, missingness

Read when computing replicate correlation, a sample-swap check, CV, or a completeness filter, or when choosing how to impute.

## Replicate Correlation on log2

**Goal:** Quantify reproducibility without letting a few abundant proteins fake agreement.

**Approach:** Correlate on log2 intensities (variance-stabilized, high-abundance tail compressed), report within-group pairs, and flag a sample correlating better with another group as a possible swap.

```python
import sys; sys.path.insert(0, "scripts")
from matrix_metrics import replicate_correlation, cross_group_correlation  # scripts/matrix_metrics.py
```

Summarize only `status == 'measured'` rows: a singleton group appears as `not_measurable_n1` with r = NaN, which means UNCHECKED, not clean. A swapped pair still correlates 0.93-0.96 with its own group's other members, so within-group r cannot show a swap; `cross_group_correlation` reports each sample's mean centred r to its own group and to every other group and sets `possible_swap` when another group matches best. A flagged sample is a swap or relabel candidate, not proof: confirm against the sample sheet.

Technical replicates r > 0.98 (instrument noise only); biological r ~ 0.90-0.98 (genuine variance, lower is expected and correct); soft floor r > 0.8 to retain a biological replicate. A Spearman check is a robustness aid only -- ranks discard the magnitude that quant QC cares about.

## Coefficient of Variation on the Linear Scale

**Goal:** Summarize per-condition precision with a number that means what it says.

**Approach:** Compute CV = SD/mean on LINEAR (non-log) intensities; if only logged values exist use the geometric-CV formula. Report the median CV per condition (the per-protein distribution is right-skewed).

```python
from matrix_metrics import median_cv_linear, geometric_cv_from_log  # scripts/matrix_metrics.py
```

Applying the base formula to log-transformed data is meaningless (Brenes 2024; see "CV computed on log-transformed data"). A group with one sample returns `status = 'not_measurable_n1'` and NaN: not measurable, not low. State normalization state, transform, and software params or the CV is uninterpretable: DIA-NN "High precision" mode silently median-normalizes, halving median CV vs "High accuracy". Technical median CV < ~10-20%, biological ~20-40%; a LOWER CV is not automatically better (loose FDR or faulty MS1 extraction produce artificially low CVs).

## Missingness Mechanism and Completeness

**Goal:** Decide how to impute by first deciding why values are missing.

**Approach:** Diagnose the missingness profile -- left-tail concentration means MNAR (left-censored, abundance-dependent), all-abundance scatter means MCAR -- and filter on completeness before imputing only the shallow remainder.

```python
from matrix_metrics import missingness_profile, completeness_filter  # scripts/matrix_metrics.py
```

kNN-imputing a genuinely-absent (MNAR) value invents mid-range abundance and KILLS a real present/absent difference; a left-shifted draw (Perseus down-shifted normal, downshift=1.8 SD below the observed mean, width=0.3 of observed SD) on an MCAR gap FABRICATES a false low and inflates a difference. Match the imputer to the mechanism. The imputation mechanics themselves are quantification.
