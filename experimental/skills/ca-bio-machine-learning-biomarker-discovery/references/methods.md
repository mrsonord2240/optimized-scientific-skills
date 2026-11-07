# Selection methods and limitations

Historical version and empirical claims are unverified. Feature selection, scaling and tuning must be nested within appropriate independent training splits for performance estimation.

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

`the original source lasso demonstration` runs the leakage and stability demonstrations (with the permuted-label and rescaling checks) on synthetic data; `the original source Boruta demonstration` compares Boruta with a cross-validated L1 fit on a co-expressed module. Both run from `scripts/`.

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
