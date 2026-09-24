---
name: bio-crispr-screens-hit-calling
category: Data Analysis
description: Cross-method decision tree for calling hits in pooled CRISPR screens. Catalogs statistical models (MAGeCK RRA, MAGeCK MLE, BAGEL2, drugZ, JACKS, Chronos, CERES), experimental designs each is built for, failure modes outside design domain, reconciliation when methods disagree, multiple-testing and effect-size thresholds, the order of operations (count -> QC -> CN-correct -> hit-call -> validate), the second-best-sgRNA conservative rule, and consensus-hit strategy. Use when choosing among MAGeCK / BAGEL2 / drugZ / JACKS / Chronos for a given design, reconciling disagreement across two or three methods on the same screen, deciding whether to require consensus, gating downstream validation by hit-confidence tier, or interpreting unstable hit lists across reruns.
tool_type: mixed
primary_tool: MAGeCK
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: MAGeCK 0.5.9+, BAGEL2 2.0, drugZ Aug 2019+, JACKS 0.2.0+, Chronos 2.0+ (DepMap), CERES 1.0+, pandas 2.2+, numpy 1.26+, scipy 1.12+, statsmodels 0.14+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mageck --version`, `BAGEL.py version`, `python drugz.py --help`
- Python: `pip show crispr_chronos` (JACKS installs from GitHub, not PyPI)

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Hit Calling Decision Tree

**"Identify significant hits in my CRISPR screen"** -> Choose the analysis method that matches the experimental design, statistical assumptions, and quality grade of the screen. Reconcile across methods when high-stakes hits must be validated.

The primary hit-calling methods cover non-overlapping niches; the decision is not "which is best" but "which matches the design." Model comparison and design rationale: [references/method-catalog.md](references/method-catalog.md).

| Design / question | Primary method | Why | Secondary check |
|--------------------|----------------|-----|------------------|
| Two-condition essentiality, one cell line, no CN concerns | MAGeCK RRA | Robust, fast, gold-standard for ranked analysis | BAGEL2 (Bayes factor on same data) |
| Time course (3+ timepoints) | MAGeCK MLE | RRA cannot model multi-condition | JACKS (efficacy-aware) |
| Multi-cell-line panel (cancer dependency) | Chronos | Models CN bias + screen quality jointly | MAGeCK MLE per line + meta-analysis (run per line first, pool downstream; a joint MLE without per-line indicator covariates dilutes per-line signal) |
| Drug screen (vehicle vs drug) | drugZ | Bidirectional Z; vehicle-anchored | MAGeCK MLE with dose covariate |
| Multi-screen joint, same library | JACKS | Shared efficacy; enables ~2.5x smaller screens | MAGeCK MLE; results should converge |
| Essentiality classification with reference sets | BAGEL2 | Bayes factor with CEGv2/NEGv1 calibration | MAGeCK RRA |
| Combinatorial / paired guide | MAGeCK MLE with GI scoring | Models interaction term; see [[combinatorial-screens]] | Custom GI scoring |
| Single-cell perturbation (Perturb-seq) | SCEPTRE | NB GLM + permutation; see [[perturb-seq-analysis]] | Mixscape pre-filter |
| Cancer-line copy-number screen | Chronos (preferred) or CERES | Joint CN-bias + gene-effect modeling; see [[copy-number-correction]] | CRISPRcleanR pre-hoc + MAGeCK |
| Cancer line + multi-batch | Chronos | Models CN and batch jointly | MAGeCK MLE with batch covariate |
| Variant function (base / prime editing) | Custom + CRISPResso2 | Editing outcomes, not guide dropout; see [[base-editing-analysis]] | -- |

## Run All Five on the Same Data (Consensus Strategy)

**Goal:** For high-stakes hits (drug-target nomination, paper-level claims), require agreement across 2-3 orthogonal methods.

**Approach:** Run MAGeCK + BAGEL2 + (drugZ or JACKS) on the same count matrix; rank by each; classify hits as called by 1, 2, or 3 methods.

```bash
python scripts/consensus_hits.py mageck.gene_summary.txt bagel.bf.txt drugz.txt -o consensus.tsv
```

`scripts/consensus_hits.py` also exposes `consensus_hits()` and `_check_comparable()` for import. It merges the three tables on gene, applies the thresholds below (FDR<0.05, BF>6, drugZ `fdr_synth`<0.05), adds `consensus_count`, and warns when a pair of hit sets is not enriched for overlap (the mismatched-comparison signature; see [references/failure-modes.md](references/failure-modes.md)). The two-method MAGeCK + BAGEL2 version is `examples/consensus_hits.py`.

**Confidence tiers:**

| Tier | Definition | Validation requirement |
|------|------------|-------------------------|
| Tier 1 (high) | Called by 3/3 methods | Arrayed validation; orthogonal modality (CRISPRi if originally Cas9) |
| Tier 2 (medium) | Called by 2/3 methods | Arrayed validation in matched line |
| Tier 3 (exploratory) | Called by 1/3 methods | Treat as hypothesis; further screens before publication |

## Reconciliation: When Two Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MAGeCK significant, BAGEL2 not | BAGEL2 trained on CEGv2/NEGv1; gene is essential but not in reference | Trust MAGeCK; flag for follow-up |
| BAGEL2 significant, MAGeCK not | BAGEL2 has tumor-suppressor sensitivity MAGeCK lacks | Investigate sgrna_summary for one weak guide |
| MAGeCK significant, JACKS not | JACKS down-weighted one outlier guide | Trust JACKS if guides agree; outlier may be off-target |
| Chronos and MAGeCK disagree on cancer line | Chronos accounts for CN; MAGeCK does not | Trust Chronos; apply [[copy-number-correction]] |
| drugZ significant, MAGeCK not on drug screen | drugZ bidirectional Z is more sensitive | Trust drugZ for chemogenomic; MAGeCK may miss small effects |
| MAGeCK MLE significant, MAGeCK RRA not in 2-condition | Beta-score effect size is significant but rank-based not | Trust MLE if guides consistent; RRA may be over-conservative |
| All methods disagree | Either no real biology or all methods are mis-applied | Stop. Re-audit QC; check chemistry / library / design matrix |
| BAGEL2 BF changes between reruns on identical input -- not a real method disagreement, but easy to mistake for one | BAGEL2's `bf` step defaults to a clock-derived random seed; two unseeded runs on the same real HAP1 TKOv3 data differed by up to 26.7 BF and flipped 33/18,053 genes across BF>6 | Always pass a fixed `-s <int>` seed to `fc`/`bf`/`pr` and verify two reruns are byte-identical before trusting any single BF table or treating a rerun difference as new biology -- see [[bagel-essentiality]]'s "Reproducibility: Fixing the Random Seed" section |

## Multiple-Testing Correction Conventions

| Method | Native correction | Cross-method comparison |
|--------|---------------------|---------------------------|
| MAGeCK RRA | BH per direction | `neg|fdr`, `pos|fdr` |
| MAGeCK MLE | BH per condition | `<cond>|fdr` |
| BAGEL2 | Bootstrap BF; reports BF threshold | BF > 6 ≈ 90% posterior (Hart 2017); ~5% FDR by convention |
| drugZ | BH per direction | `fdr_synth`, `fdr_supp` |
| JACKS | Posterior probability + BH | `fdr_log10` (log10 FDR) |
| Chronos | DepMap gene-effect probability | `effect_probability` |

**Reconciliation:** BF >6 in BAGEL2 corresponds to ~90% posterior probability (Hart 2017 G3, by overlap with CEGv2) and is commonly used as a stringent cutoff roughly comparable to MAGeCK FDR 0.05. Treat that equivalence as an approximate convention, not an exact calibration. drugZ FDR is per-direction; the `fdr_synth` and `fdr_supp` columns are independent BH corrections.

## Correlating MAGeCK and BAGEL2 Scores (Sign/Scale Caution)

`neg|score` (MAGeCK RRA) is p-value-like: smaller = more essential. `BF` (BAGEL2) is a
log-likelihood ratio: larger = more essential. The two statistics run in opposite directions on
the same biology. Computing Spearman rho directly between the raw columns therefore produces a
strongly *negative* number even when the methods agree strongly -- sign-correct one column first
(use `-neg|score`) before comparing rho to a disagreement threshold.

**Worked example, real HAP1 TKOv3 T0-vs-T18 data (n=18,053 genes with both scores):**

| Comparison | Spearman rho |
|---|---|
| `neg|score` vs `BF`, uncorrected | -0.806 |
| `-neg|score` vs `BF`, sign-corrected | +0.806 |

Both rows describe the same real agreement between MAGeCK and BAGEL2 on this screen -- only the
sign flips. Applying a literal "rho <0.6 means disagreement" rule to the uncorrected row would
falsely flag two methods that in fact concur strongly. Always sign-correct before comparing rho
to a threshold, and note it in whatever you report ("rho = 0.81 after sign-correcting `neg|score`
to make both statistics increase with essentiality").

## Order of Operations

```
1. Library design (see library-design)         <- design quality dictates hit calling
2. Plasmid pool sequencing                     <- baseline; non-negotiable
3. Run screen at MOI 0.3, 500x coverage
4. Sequence endpoint
5. Run mageck count                            <- generates raw + normalized counts
6. Screen QC (see screen-qc)                   <- gates downstream method choice; run the CEGv2/NEGv1 PR-AUC first
                                                  (PR-AUC <0.5 = no signal however many hits MAGeCK calls; >0.7 to interpret)
7. Copy-number correction if cancer line       <- CRISPRcleanR or Chronos; see copy-number-correction
8. Batch correction if multi-batch             <- see batch-correction; batch covariates in MAGeCK MLE, or Chronos
9. Hit calling (this skill)                    <- choose method by design
10. Consensus across 2-3 methods               <- for high-stakes hits
11. Orthogonal validation                      <- arrayed; different chemistry
12. Pathway analysis                           <- see pathway-analysis/gsea
```

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| MAGeCK RRA FDR (gene-level) | <0.05 | Li 2014; standard publication |
| MAGeCK RRA LFC | abs(LFC) >1 | 2-fold; biological |
| BAGEL2 Bayes Factor | >6 standard; >12 stricter | Hart 2017; BAGEL convention |
| drugZ FDR | <0.05 per direction | Colic et al. 2019 |
| JACKS fdr_log10 | <-1 (FDR <0.1); <-2 (FDR <0.01) | Standard FDR convention |
| Chronos dependency probability | >0.5 | DepMap convention (dependency-probability cutoff) |
| Tier 1 consensus (3 methods) | 100% agreement | High confidence; minimal validation needed |
| Tier 2 consensus (2 of 3) | 67% agreement | Arrayed validation required |
| Tier 3 (1 method only) | Hypothesis; flag for follow-up | Multiple screens or arrayed required |
| Second-best sgRNA rule | Second-best LFC also passes threshold | Reduces single-guide outliers |

This table is the canonical source for MAGeCK-FDR/BAGEL-BF defaults: `scripts/consensus_hits.py` and
`examples/consensus_hits.py` must both use FDR<0.05/BF>6 to match it. On real HAP1
TKOv3 data, FDR<0.05/BF>6 gives 844 Tier-1 consensus genes; the looser FDR<0.1/BF>5 pairing gives
1131 (+34%) -- pick one number, not whichever default a given script happens to hardcode.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| All genes significant in MAGeCK RRA | Heavy selection breaks median norm | `--norm-method control`; or use BAGEL2 |
| BAGEL2 returns no hits despite known essentials | Wrong reference gene set | Verify CEGv2/NEGv1 files match library |
| drugZ output empty | Used Day 0 as control instead of vehicle | Re-run with vehicle as control |
| Chronos errors out | Missing CN profile for cell line | Use CRISPRcleanR (unsupervised) instead |
| Methods disagree by orders of magnitude | Quality issue or design mismatch | Re-audit QC; reconcile via tier consensus |
| Empty tier 1 consensus | No real biology, QC failure, OR merged files are not from the same experimental comparison | Check per-method QC (screen-qc) and `_check_comparable()`'s overlap-enrichment warning first; verify all inputs are the same comparison; only then re-audit QC |
| Single-guide-driven hits, or genes with only 1 sgRNA in the library | Outlier sgRNA, or library design has no second guide to check | Apply second-best rule ([references/second-best-and-custom-zscore.md](references/second-best-and-custom-zscore.md)); treat `single_guide=True` the same as a failed check; orthogonal validate |
| BAGEL2 BF differs between two runs on the same input | `bf` step unseeded by default (clock-derived) | Always pass `-s <fixed-int>`; see [[bagel-essentiality]] Reproducibility section |

## Reference Files

| File | Read when |
|------|-----------|
| [references/method-catalog.md](references/method-catalog.md) | Comparing the seven methods' models, RRA vs MLE, or why each was built |
| [references/second-best-and-custom-zscore.md](references/second-best-and-custom-zscore.md) | Filtering hits for single-guide outliers, or calling hits with a custom z-score |
| [references/failure-modes.md](references/failure-modes.md) | Methods disagree by 200+ hits, or the consensus is empty |

## References

- Li W et al. 2014. *Genome Biol* 15:554. MAGeCK alpha-RRA.
- Li W et al. 2015. *Genome Biol* 16:281. MAGeCK MLE.
- Kim E & Hart T. 2021. *Genome Med* 13:2. BAGEL2.
- Colic M et al. 2019. *Genome Med* 11:52. drugZ.
- Allen F et al. 2019. *Genome Res* 29:464. JACKS.
- Dempster J et al. 2021. *Genome Biol* 22:343. Chronos.
- Meyers R et al. 2017. *Nat Genet* 49:1779. CERES.
- Hart T & Moffat J. 2016. *BMC Bioinformatics* 17:164. BAGEL Bayes factor framework.
- Hart T et al. 2017. *G3* 7:2719. CEGv2/NEGv1 calibration.

## Related Skills

- crispr-screens/mageck-analysis - Full MAGeCK RRA + MLE detail
- crispr-screens/bagel-essentiality - Full BAGEL2 detail
- crispr-screens/drugz-chemogenomic - Full drugZ detail for drug screens
- crispr-screens/jacks-analysis - Full JACKS detail and library calibration
- crispr-screens/copy-number-correction - Chronos, CERES, CRISPRcleanR
- crispr-screens/screen-qc - Quality gates that drive method choice
- crispr-screens/library-design - Library type dictates analysis method
- crispr-screens/combinatorial-screens - GI scoring (synthetic lethality)
- crispr-screens/perturb-seq-analysis - SCEPTRE for single-cell screens
- crispr-screens/batch-correction - Multi-batch normalization upstream of hit calling
- pathway-analysis/gsea - Downstream pathway enrichment
- pathway-analysis/go-enrichment - GO enrichment of hit lists
