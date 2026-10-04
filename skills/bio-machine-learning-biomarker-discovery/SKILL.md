---
name: bio-machine-learning-biomarker-discovery
category: Data Analysis
description: Use when identifying candidate biomarkers from high-dimensional omics data, deciding between an all-relevant and a minimal-optimal selector, or judging whether a selected gene set is reproducible.
tool_type: python
primary_tool: boruta
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: numpy 2.5, pandas 3.0, scikit-learn 1.9, boruta 0.4.3, mrmr-selection 0.2.8 (the bundled scripts were run against these).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

BorutaPy expects numpy arrays (pass `X.values`); boruta 0.4+ works on current numpy, while 0.3.x breaks on the removed `np.float`/`np.int` aliases. scikit-learn 1.8 deprecated the `penalty=` argument of `LogisticRegression`/`LogisticRegressionCV` (removed in 1.10): leave `penalty` unset and express the penalty with `l1_ratio` (`l1_ratios` in the CV class) and `C`; `l1_ratio=1` is L1, `0` is L2, in between is elastic net. If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Feature Selection for Biomarker Discovery

**"Find the biomarkers in my omics data"** -> First decide which question is being answered (all-relevant vs minimal-optimal), then select features INSIDE a resampling loop, then quantify stability -- because a selected list means little without it.
- All-relevant (which genes carry signal?): `BorutaPy(rf)`
- Minimal-optimal (smallest predictive set?): `ElasticNetCV`, `LogisticRegressionCV(solver='saga', l1_ratios=[...])`
- Stability (does the list reproduce?): bootstrap selection frequencies + a stability index

## The Single Most Important Modern Insight -- Most Gene Signatures Do Not Replicate, and Significance Is the Wrong Bar

A signature being "significantly associated with outcome" is near-worthless evidence: *random* gene sets -- and signatures of biologically irrelevant phenomena -- are significantly associated with breast-cancer survival, often matching published prognostic signatures, because the transcriptome is dominated by a few axes (proliferation) that almost any large gene set captures (Venet 2011). The correct null is not "no association" but *random gene sets of equal size* plus a proliferation meta-gene. Two further hard facts complete the picture: many disjoint gene lists predict equally well (Ein-Dor 2005), so non-overlap with a prior list is the *expected* result, not a contradiction; and obtaining a *stable* list (as opposed to an accurate predictor) needs on the order of thousands of samples (Ein-Dor 2006), far more than typical omics n.

The operational consequences run through every section below: report a stability index next to accuracy; benchmark against a random-signature and proliferation-meta-gene null; never interpret the specific genes a minimal-optimal selector kept as "the biomarkers"; and keep selection inside the cross-validation loop or the reported performance is fiction.

## All-Relevant vs Minimal-Optimal (the distinction usually conflated)

This axis matters more than filter/wrapper/embedded. Choosing the wrong one is the most common conceptual error in applied biomarker papers.

- **Minimal-optimal** = the *smallest* subset giving optimal prediction (LASSO, RFE, forward selection). If two genes are correlated and both informative, it keeps **one and drops the other**; the dropped gene is still biologically relevant. Minimal-optimal sets are non-unique, unstable, and systematically exclude redundant-but-real features. *Absence from a minimal-optimal set is not evidence of irrelevance.*
- **All-relevant** = *every* feature carrying information, redundant or not (Boruta: keep anything beating the best "shadow" permuted feature). This is the right framing for *biological interpretation* -- the whole co-expression module is wanted, not one representative.

Decision rule: parsimonious assay with few measurements -> minimal-optimal; understand biology / enumerate implicated genes / pathway analysis -> all-relevant; stable deployable signature -> elastic net or stability selection.

## Methods Taxonomy

| Family | Method | Optimizes | Redundancy handling | Output | Key trap |
|--------|--------|-----------|---------------------|--------|----------|
| Filter (univariate) | t-test / `SelectKBest(f_classif)` | Marginal association, one gene at a time | None (keeps correlated blocks) | Ranked list | Ignores multivariate structure; huge multiplicity |
| Filter (multivariate) | mRMR (Peng 2005) | Relevance minus redundancy | Explicit penalty | Ranked K | Greedy/first-order; K must still be chosen |
| Wrapper | RFE / RFECV; SVM-RFE | A specific model's accuracy | Indirect | Ranked subset | Expensive; **must be inside CV**; SVM-RFE needs a linear kernel |
| Embedded | LASSO (Tibshirani 1996) | Prediction + L1 sparsity | **None** -- arbitrarily keeps one of a correlated group | Sparse coefs | Unstable under collinearity; caps at n features when p>n |
| Embedded | Elastic net (Zou-Hastie 2005) | Prediction + L1+L2 grouping | Keeps correlated groups together | Sparse coefs | Two hyperparameters; still not "causal" |
| All-relevant | Boruta (Kursa 2010) | Every feature beating shadow features | Keeps all relevant (redundant included) | Confirmed/Tentative/Rejected | Slow; returns redundant sets by design |
| Meta / stability | Stability selection (Meinshausen 2010; Shah-Samworth 2013) | Selection probability under subsampling | Inherits base learner | Selection frequencies + threshold | Error bounds assume exchangeability omics violates |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Want every implicated gene for pathway/biology interpretation | Boruta (all-relevant), or stability-based consensus | Keeps whole correlated modules, not one representative |
| Want a small deployable assay/signature | Elastic-net (not bare LASSO); report stability | L2 grouping keeps correlated genes together and resamples more stably |
| p is huge (>20k); selection is slow | Univariate pre-filter to a few thousand, then Boruta/elastic-net, all inside the CV fold | Cheap dimensionality cut; never pre-filter on the full dataset |
| Need to report model performance | Wrap selection in a `Pipeline`, estimate by nested CV | Selection outside CV inflates AUC to ~perfect on pure noise |
| Single-cell biomarker across conditions | Pseudobulk per donor, then select at the donor level | The unit is the donor, not the cell (Squair 2021); cells are pseudoreplicates |
| Want to know which genes "drive" a trained model | -> machine-learning/prediction-explanation | SHAP ranking is not validated selection |
| Want unbiased accuracy/calibration of the selected model | -> machine-learning/model-validation | Selection is one step; validation is its own discipline |

## Leakage-Safe Selection (the single most damaging error to avoid)

**Goal:** Estimate the performance of a selection-plus-model pipeline without optimistic bias.

**Approach:** Put selection in a `Pipeline` so it is re-fit on each training fold only; the held-out fold never informs which features are kept. Selecting the top-k features on the *whole* dataset before cross-validating the classifier produces near-zero apparent error even on pure noise (Ambroise-McLachlan 2002). Selection is where almost all overfitting capacity lives when p>>n.

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold

pipe = Pipeline([
    ('scale', StandardScaler()),                            # fit on training folds only
    ('select', SelectKBest(f_classif, k=20)),               # re-fit per fold -> no leakage
    ('clf', LogisticRegression(max_iter=5000)),
])
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=0)
auc = cross_val_score(pipe, X, y, cv=cv, scoring='roc_auc')   # held-out folds, k fixed in advance
print(f'In-pipeline CV AUC (k fixed): {auc.mean():.3f} +/- {auc.std():.3f}')
```

This is single-level CV: valid only because `k` was not tuned on these folds. If `k` (or any hyperparameter) is tuned, nest the search so the outer folds never see the tuning:

```python
from sklearn.model_selection import GridSearchCV

tuned = GridSearchCV(pipe, {'select__k': [10, 20, 50]}, scoring='roc_auc',
                     cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=1))
nested_auc = cross_val_score(tuned, X, y, cv=cv, scoring='roc_auc')
```

Univariate selection is not the only in-fold option. mRMR needs a DataFrame/Series, so wrap it as a transformer whose `fit` runs on the training fold only:

```python
import numpy as np
import pandas as pd
from mrmr import mrmr_classif
from sklearn.base import BaseEstimator, TransformerMixin

class MRMRSelect(BaseEstimator, TransformerMixin):
    def __init__(self, K=20):
        self.K = K
    def fit(self, X, y):
        self.cols_ = mrmr_classif(X=pd.DataFrame(X), y=pd.Series(np.asarray(y)), K=self.K, show_progress=False)
        return self
    def transform(self, X):
        return np.asarray(X)[:, self.cols_]

mrmr_pipe = Pipeline([('scale', StandardScaler()), ('select', MRMRSelect(K=20)),
                      ('clf', LogisticRegression(max_iter=5000))])
mrmr_auc = cross_val_score(mrmr_pipe, X, y, cv=cv, scoring='roc_auc')
```

The standalone Boruta/elastic-net blocks below select features on a full matrix to *discover* candidates; that is fine for discovery, but any performance number must come from a Pipeline pattern above, with scaling and selection inside the fold.

## All-Relevant: Boruta

**Goal:** Enumerate every feature carrying signal, including redundant co-expressed genes.

**Approach:** Compare each real feature's importance to the maximum importance of permuted "shadow" features over many iterations; confirm features that consistently beat the best shadow.

```python
from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier

rf = RandomForestClassifier(n_estimators=100, n_jobs=-1, class_weight='balanced', max_depth=5, random_state=42)
# perc=100 uses the max shadow importance (strict); two_step (default True) controls the multiple-testing correction.
boruta = BorutaPy(rf, n_estimators='auto', perc=100, two_step=True, max_iter=100, random_state=42)
boruta.fit(X.values, y.values)                              # numpy arrays, not pandas

confirmed = X.columns[boruta.support_]                      # all-relevant set (redundant by design)
tentative = X.columns[boruta.support_weak_]
```

## Minimal-Optimal: Elastic Net (prefer over bare LASSO)

**Goal:** A small, stable predictive signature from correlated omics features.

**Approach:** Use elastic net, whose L2 term induces a grouping effect so correlated genes enter or leave together; standardize first because the penalty is scale-sensitive. Bare LASSO keeps one arbitrary member of a correlated group and flips on tiny data perturbations.

```python
from sklearn.linear_model import LogisticRegressionCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# Scaler and penalty search live in one pipeline, so the scaler sees training rows only.
# saga is the only solver supporting mixed L1/L2 (0 < l1_ratio < 1); C = 1/lambda (opposite of alpha in Lasso/ElasticNet).
enet = make_pipeline(
    StandardScaler(),
    LogisticRegressionCV(solver='saga', l1_ratios=[0.1, 0.5, 0.9], Cs=20, cv=5, max_iter=10000,
                         scoring='accuracy', use_legacy_attributes=False))
enet.fit(X, y)
selected = X.columns[enet[-1].coef_[0] != 0]
```

**Scoring sets the signature size, so it is pinned.** `LogisticRegressionCV` scores its inner folds with accuracy unless told otherwise, and scikit-learn plans to change that default (1.11); leave `scoring` explicit. On Golub ALL/AML (72 samples, 1,000 highest-variance probes, one 5-fold partition) the selected-probe count was 601 for `neg_log_loss`, 203 for `roc_auc` and 148 for `accuracy`, and the held-out AUC was within 0.01 across all three (0.98-0.99, Golub is near-separable). With labels permuted, `neg_log_loss` still selected 16 probes while the other two selected none. For a compact signature the default here is `accuracy`; `neg_log_loss` favours probability fit and keeps many more features. Measured only on this one dataset and partition: re-check the size on your own data, and treat accuracy as coarse on small or imbalanced samples.

**Cost.** Measured on Golub (72 samples, 24-thread Windows machine, another job running): this elastic-net block took about 34 s on the 1,000 probes and about 263 s on all 7,129; pre-filter before running it on more.

## Stability: Are the Selected Features Reproducible?

**Goal:** Distinguish a robust signature from a resampling accident, and report stability alongside accuracy.

**Approach:** Run an L1 selector on many n/2 subsamples, count per-feature selection frequency, keep features above a threshold (0.6 is the common default), and compute a chance-corrected stability index (Nogueira 2018: handles variable-size selections; the older Kuncheva index needs equal-size subsets and breaks for LASSO). Two choices make the result trustworthy rather than unit-dependent:

- **Scale inside every resample, and fix the budget instead of C.** A fixed `C` means a different number of features for every choice of units (on raw Golub probes `C=0.1` kept ~100 per subsample and reported dozens of "stable" probes even for permuted labels). `scripts/stability_selection.py` standardizes each subsample, walks the L1 path from the strongest penalty, and stops at `q` selected features per subsample (Meinshausen-Buhlmann 2010). Under exchangeability the expected number of false selections is at most `q^2 / ((2*pi_thr - 1) * p)`; set the tolerated expected count `EV` (1 below) and solve for `q = sqrt(EV * (2*pi_thr - 1) * p)`.
- **Report the permuted-label null next to the real count.** Re-run the identical procedure on shuffled labels; the null count should be near zero, and a real count that does not clearly exceed it is not a signature. Stability is not evidence of signal: the index can be high for a feature set that is stably wrong.

The function is deterministic: the same data, `q` and `seed` give identical frequencies (checked twice on raw and rescaled input, 1,000 and 7,129 probes). Unit invariance is approximate, because rescaling perturbs the path's floating-point arithmetic and a marginal feature can flip; it held in the runs below. On Golub ALL/AML with `pi_thr` 0.6, `EV` 1 and `seed=0`:

| Probes (`q`) | Stable, real labels (raw / standardized / rescaled units) | Stable, permuted labels |
|---|---|---|
| 1,000 highest-variance (14) | 2 / 2 / 2 | 0 in 10 of 11 permutations, 1 stable probe (frequency 0.81) in one |
| all 7,129 (37) | 5 / 5 / 5 | 0 in all 11 permutations (highest frequency 0.57) |

A non-zero null count is possible, so run several permutations and report the count over runs, not a single 0. One call (100 subsamples) took 1-6 s at 1,000 probes and 15-24 s at 7,129, so the repeats add up.

```python
import sys
import numpy as np
sys.path.insert(0, 'scripts')                              # run from the Skill directory
from stability_selection import stability_selection        # deterministic for a given seed; returns (frequencies, Nogueira index)

pi_thr, EV = 0.6, 1                                        # pi_thr: Meinshausen-Buhlmann default; EV: tolerated false selections
q = int(np.sqrt(EV * (2 * pi_thr - 1) * X.shape[1]))       # per-subsample budget
freq, stability = stability_selection(X, y, q=q, seed=0)
stable = X.columns[freq > pi_thr]
freq_null, _ = stability_selection(X, np.random.default_rng(1).permutation(y), q=q, seed=0)
print(f'{len(stable)} stable features on real labels vs {(freq_null > pi_thr).sum()} on permuted labels')
print('no features selected' if np.isnan(stability) else f'Nogueira stability = {stability:.2f}')
```

`python scripts/lasso_biomarker.py` runs the leakage and stability demonstrations (with the permuted-label and rescaling checks) on synthetic data; `python scripts/boruta_feature_selection.py` compares Boruta with a cross-validated L1 fit on a co-expressed module. Both run from `scripts/`.

## Failure Modes

Read [`references/failure-modes.md`](references/failure-modes.md) when a selected list looks unstable, a signature fails to replicate, or selectors disagree.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Samples for a *stable* gene list ~ thousands | Ein-Dor 2006 | Small effects need large n for reproducible membership (accuracy needs far fewer) |
| Selection inside every CV fold; nested CV for tuning | Ambroise 2002; Simon 2003 | Selection outside CV gives ~0% error on noise |
| Stability threshold pi_thr ~ 0.6-0.9 | Meinshausen-Buhlmann 2010 | Selection-frequency cutoff; tune to false-positive cost |
| Random-signature null | Venet 2011 | Benchmark against size-matched random sets + proliferation meta-gene |
| Single-cell unit = donor (pseudobulk) | Squair 2021 | Cells are pseudoreplicates |
| Biomarker clinical translation rate <1% | Kern 2012 | Sets expectations; failures follow a foreseeable taxonomy |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| BorutaPy raises on pandas input (or on `np.float` with boruta 0.3.x) | Needs numpy arrays; 0.3.x predates numpy's alias removal | Pass `X.values`, `y.values`; use boruta 0.4+ |
| Regularization strength backwards | `C=1/lambda` (logistic) vs `alpha` (Lasso/ElasticNet) are opposite conventions | Verify which API; small C = strong shrinkage |
| Mixed-penalty solver error, or `penalty=` FutureWarning (sklearn 1.8+) | Only `solver='saga'` supports `0 < l1_ratio < 1`; `penalty=` is deprecated | Set `solver='saga'`, pass `l1_ratio(s)`, leave `penalty` unset |
| `mrmr_classif` returns wrong type | Pandas backend needs a DataFrame X and Series y | Pass `X` DataFrame, `y=pd.Series(y)`; K must still be chosen |
| glmnet signature unstable across runs (R) | Used `lambda.min` | Use `lambda.1se` for a sparser, more reproducible set |

## References

- Tibshirani R. 1996. Regression shrinkage and selection via the lasso. *J R Stat Soc B* 58:267-288.
- Goring HHH, Terwilliger JD, Blangero J. 2001. Large upward bias in estimation of locus-specific effects from genomewide scans. *Am J Hum Genet* 69:1357-1369.
- Ambroise C, McLachlan GJ. 2002. Selection bias in gene extraction on the basis of microarray gene-expression data. *PNAS* 99:6562-6566.
- Simon R, Radmacher MD, Dobbin K, McShane LM. 2003. Pitfalls in the use of DNA microarray data for diagnostic and prognostic classification. *J Natl Cancer Inst* 95:14-18.
- Ein-Dor L, Kela I, Getz G, Givol D, Domany E. 2005. Outcome signature genes in breast cancer: is there a unique set? *Bioinformatics* 21:171-178.
- Peng H, Long F, Ding C. 2005. Feature selection based on mutual information: criteria of max-dependency, max-relevance, and min-redundancy. *IEEE Trans Pattern Anal Mach Intell* 27:1226-1238.
- Zou H, Hastie T. 2005. Regularization and variable selection via the elastic net. *J R Stat Soc B* 67:301-320.
- Ein-Dor L, Zuk O, Domany E. 2006. Thousands of samples are needed to generate a robust gene list for predicting outcome in cancer. *PNAS* 103:5923-5928.
- Kursa MB, Rudnicki WR. 2010. Feature selection with the Boruta package. *J Stat Softw* 36:1-13.
- Meinshausen N, Buhlmann P. 2010. Stability selection. *J R Stat Soc B* 72:417-473.
- Venet D, Dumont JE, Detours V. 2011. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput Biol* 7:e1002240.
- Kern SE. 2012. Why your new cancer biomarker may never work: recurrent patterns and remarkable diversity in biomarker failures. *Cancer Res* 72:6097-6101.
- Shah RD, Samworth RJ. 2013. Variable selection with error control: another look at stability selection. *J R Stat Soc B* 75:55-80.
- Nogueira S, Sechidis K, Brown G. 2018. On the stability of feature selection algorithms. *J Mach Learn Res* 18:1-54.
- Squair JW, Gautier M, Kathe C, et al. 2021. Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692.

## Related Skills

- machine-learning/model-validation - Nested CV and leakage-safe estimation of the selected model
- machine-learning/prediction-explanation - Why SHAP rankings are not a validated selection method
- machine-learning/omics-classifiers - Build a classifier from the selected features
- differential-expression/de-results - Pre-filter candidates with differential expression
- experimental-design/multiple-testing - FDR control and why it is orthogonal to selection stability
- experimental-design/power-analysis - Sample size for a stable signature vs an accurate predictor
- pathway-analysis/go-enrichment - Functional enrichment of an all-relevant gene set
