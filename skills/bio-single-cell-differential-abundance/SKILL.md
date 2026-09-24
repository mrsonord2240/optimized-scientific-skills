---
name: bio-single-cell-differential-abundance
description: Test whether cell-type proportions or composition changed between conditions in single-cell data using Milo (miloR), scCODA, sccomp, and propeller. Use when comparing cell-type proportions / composition between conditions, asking which populations expanded or contracted with treatment or disease, running neighborhood-level (cluster-free) abundance testing, or guarding against compositional shifts that masquerade as differential expression.
tool_type: mixed
primary_tool: Milo
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: miloR 2.0+, scCODA 0.1.9+, sccomp 1.8+, speckle 1.0+
Checked by execution on miloR 2.2.0, speckle 1.6.0 and scCODA 0.1.9 (R 4.4.3); the sccomp block was checked against the upstream README only, not run.

```r
BiocManager::install(c('miloR', 'sccomp', 'speckle'))   # sccomp also needs CmdStan (cmdstanr::install_cmdstan())
```

```bash
pip install scanpy sccoda   # scCODA 0.1.9 needs arviz<1 (arviz 1.x removed the arviz.data API it calls) and tf_keras
```

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Differential Abundance Testing

**"Did cell-type proportions change between conditions?"** -> Test whether populations expanded or contracted between groups, accounting for the fact that proportions are not independent.
- R (cluster-free): `miloR` - build kNN graph, define neighborhoods, `testNhoods()` with a GLM and SpatialFDR
- Python/R (cluster-based): `scCODA` (Bayesian Dirichlet-multinomial), `sccomp` (Bayesian, outlier-robust), `propeller` (speckle, arcsin-sqrt + limma)

## Governing principle

Composition data live on a SIMPLEX: proportions sum to 1, so they are NOT independent - when one population expands, every other proportion is mechanically forced down even if its absolute count never changed. Running a per-cluster t-test (or Wilcoxon) on proportions across samples is therefore invalid: it ignores the negative correlation the constraint imposes, treats each cell type as a free measurement, and produces correlated false positives (one true expansion drags down the rest, which then test as spurious "depletions"). Valid methods model the joint composition: either a Dirichlet-multinomial / log-ratio model with a reference (scCODA, sccomp) or a variance-stabilizing transform plus a linear model (propeller), or they sidestep hard clusters entirely by testing abundance on the kNN graph (Milo).
Replicates are samples, not cells. The unit of replication for a composition claim is the biological sample/donor; thousands of cells from one donor are one draw. Differential abundance needs biological replicates per condition (Milo, scCODA, sccomp, propeller all model sample-level counts), and few replicates (n<3-4/group) leave abundance shifts underpowered and unstable - more donors help, more cells per donor barely do. With n=1 per condition the donor is perfectly confounded with condition: the effect is unidentifiable, not merely underpowered. propeller refuses (`No finite residual standard deviations`); scCODA re-run on two n=1 pairs (S1 vs S5 with two seeds, S2 vs S8) returned no credible effects even where NK had doubled, which is an uninformative answer, not evidence of no change; sccomp at n=1 is untested. Require >=2 (ideally 3-4) biological replicates per group before believing any abundance call.
Differential abundance and differential expression are different questions and confound each other. A pseudobulk or cluster-level "DE" signal between conditions can be pure composition: if a cluster mixes substates and treatment shifts their ratio, the aggregated profile changes although no gene changed expression in any cell - differential abundance masquerading as differential expression, invisible if only DE is run. Always pair a condition-DE analysis (single-cell/markers-annotation, differential-expression/deseq2-basics) with a differential-abundance test and interpret them jointly.
Milo's neighborhood sampling and scCODA's HMC chain are stochastic. Set a seed (`set.seed()` before `makeNhoods`; `tf.random.set_seed()` before `CompositionalAnalysis()`, not just before `sample_hmc`) and report it with the k/`prop` or HMC settings that produced the result. Verified: the same seed reproduced Milo's neighborhood set and scCODA's inclusion probabilities exactly; unseeded runs differed (Milo 272-293 neighborhoods; scCODA inclusion probabilities of unchanged types moved by up to 0.07-0.14).

## Choosing a differential-abundance method

| Method | Model | Granularity | Use when | Fails when |
|--------|-------|-------------|----------|------------|
| Milo (miloR) | NB-GLM on kNN-neighborhood counts, SpatialFDR | Cluster-free neighborhoods | Continuous/transitional states; shifts that discrete clusters hide; want sub-cluster resolution; many cells per sample (see the power note in the Milo section) | Under-powered at ~700-850 cells/sample (0 of 284 neighborhoods significant for a 2.25x NK expansion that propeller called at FDR 9e-5); results sensitive to k and `prop`; needs an integrated embedding |
| scCODA | Bayesian Dirichlet-multinomial, log-linear, reference cell type | Discrete clusters | Cluster-level testing with the simplex bias handled; want credible effects / FDR | Reference cell type mis-chosen; very few samples; HMC tuning |
| sccomp | Bayesian beta-binomial mixed model, outlier-robust | Discrete clusters | Outliers/over-dispersion present; want joint mean + variability, random effects | Small data with weak priors; longer runtime |
| propeller (speckle) | logit transform (speckle 1.6.0 default; `transform = 'asin'` for arcsin-sqrt) + limma moderated test | Discrete clusters | Fast frequentist test, several samples/group, Seurat/SCE input | Very small sample counts; ignores some compositional coupling vs Bayesian models |
| Simple proportion t-test / chi-square | Per-cluster test on proportions | Discrete clusters | Never recommended as the primary test | Always - ignores the simplex; correlated false positives |

scCODA and sccomp are cluster-based and Bayesian and report credible/FDR-controlled effects; Milo is cluster-free and catches shifts within a cell type that clustering averages away; propeller is the fast frequentist option. Lead with a cluster-based method (propeller, scCODA or sccomp) as the primary test; run Milo as the sub-cluster complement when there are enough cells per sample, and reconcile. A Milo null does not overrule a significant cluster-level call. When methods compete, verify current best practice against installed docs.

## The reference-cell-type choice in scCODA

Compositional analysis is always relative to something. scCODA fixes one cell type as the reference assumed unchanged by the covariates, and reports every other type's change relative to it; the verdict can flip with a different reference. Choose a cell type that is biologically stable and abundant across all samples, or use `reference_cell_type='automatic'` (scCODA picks a type with low dispersion present in all samples). A reference that actually changes will bias all other calls. sccomp avoids a hard reference by modeling all groups jointly; Milo avoids it via the graph.

## Adjusting for nuisance covariates and confounded designs

**Goal:** Adjust the abundance model for technical or biological nuisances (sequencing batch, timing, sex, age) and recognize when adjustment cannot help.

**Approach:** Add the nuisance as an extra additive term in the model formula with the condition of interest last; the test then reports the condition effect holding the nuisance constant. The nuisance column must vary within each condition - if a batch is perfectly confounded with condition (e.g. all controls sequenced in batch 1, all treated in batch 2), the term is unidentifiable and the test is invalid; the fix is experimental (multiplex conditions across batches), not statistical.

```r
# Milo: batch added before condition; batch column lives in design.df
design <- distinct(as.data.frame(colData(milo))[, c('sample', 'batch', 'condition')])
rownames(design) <- design$sample
da <- testNhoods(milo, design = ~ batch + condition, design.df = design, reduced.dim = 'PCA')
```

```python
# scCODA: additive patsy formula; covariate columns must be in the count table
data = dat.from_pandas(counts, covariate_columns=['sample', 'batch', 'condition'])
model = mod.CompositionalAnalysis(data, formula='batch + condition', reference_cell_type='automatic')
```

```r
# sccomp: nuisance added to formula_composition (and optionally formula_variability)
res <- sccomp_estimate(counts_tbl, formula_composition = ~ batch + condition, .sample = sample, .cell_group = cell_type, .count = count, cores = 1)
```

Diagnose confounding before modeling: cross-tabulate batch x condition; if a batch maps to a single condition, no covariate term recovers the effect. Build Milo's kNN graph on a batch-corrected embedding, but keep batch in the GLM design as well, since integration and design adjustment address different residual structure. Check that the adjusted Milo call agrees with a cluster-level test: in the audit (batch balanced across conditions, 8 samples) `~ batch + condition` turned Milo's 0 significant neighborhoods into 16, 15 of them CD14+ monocytes, a population that did not change, and only 1 NK; three parameters on 8 samples leave 5 residual df, so report an adjusted-only neighborhood call as unconfirmed.

## Milo - cluster-free neighborhood abundance (R)

**Goal:** Test differential abundance on kNN neighborhoods so shifts within and between cell types are both visible.

**Approach:** Build the Milo object from an integrated reduced dimension, sample representative neighborhoods, count cells per sample per neighborhood, then fit a GLM with `testNhoods` and control the graph-aware SpatialFDR; annotate neighborhoods back to cell types for interpretation.

```r
library(miloR)
library(SingleCellExperiment)
library(dplyr)

set.seed(42)
milo <- Milo(sce)
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = 'PCA')
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE, reduced_dims = 'PCA')
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), samples = 'sample')

design <- data.frame(colData(milo))[, c('sample', 'condition')]
design <- distinct(design)
rownames(design) <- design$sample
milo <- calcNhoodDistance(milo, d = 30, reduced.dim = 'PCA')

da <- testNhoods(milo, design = ~ condition, design.df = design, reduced.dim = 'PCA')
da <- annotateNhoods(milo, da, coldata_col = 'cell_type')
table(da$SpatialFDR < 0.1, da$cell_type)
```

**Power.** Milo's NB-GLM on ~100-cell neighborhoods is under-powered at ordinary experiment sizes. Audit run (8 donors, 4 per group, 674-837 cells per sample, NK 6% -> 13%, a 2.25x expansion in every donor): 0 of 284 neighborhoods at SpatialFDR < 0.1 and no raw p < 0.05, while propeller gave FDR 9e-5 and scCODA inclusion probability 1.0 on the same cells; widening k to 120 (median 353 cells/neighborhood) only reached a best NK raw p of 0.067, with the logFC sign correct (+1.3 to +2.3) throughout. With hundreds of cells per sample, or a population of ~5-10% of cells, make the cluster-based result the primary test and report Milo's logFC as supporting effect size, not as a null. Not tested: how many cells per sample Milo needs.

`k` and `prop` trade resolution against power: larger neighborhoods are better powered but blur fine shifts. SpatialFDR (not raw p) corrects for overlapping neighborhoods - report it. A neighborhood with a mixed `cell_type` fraction is a genuinely transitional region, not a labeling error.

## scCODA - Bayesian cluster-level composition (Python)

**Goal:** Test cluster proportion changes while handling the simplex's negative-correlation bias.

**Approach:** Build a per-sample cell-type count table with covariates, fit the Dirichlet-multinomial model against a reference cell type, sample the posterior, then read credible effects at a chosen FDR.

```python
import pandas as pd
import tensorflow as tf
from sccoda.util import cell_composition_data as dat
from sccoda.util import comp_ana as mod

counts = pd.crosstab(adata.obs['sample'], adata.obs['cell_type']).reset_index()
meta = adata.obs[['sample', 'condition']].drop_duplicates()
counts = counts.merge(meta, on='sample')

data = dat.from_pandas(counts, covariate_columns=['sample', 'condition'])
tf.random.set_seed(42)   # import tensorflow as tf. sample_hmc has no seed argument, and the initial state is drawn when the model is built, so seed BEFORE CompositionalAnalysis()
model = mod.CompositionalAnalysis(data, formula='condition', reference_cell_type='automatic')
result = model.sample_hmc()   # num_results=20000, num_burnin=5000 by default (n_burnin is not an argument; num_results must exceed num_burnin)
result.set_fdr(est_fdr=0.1)
result.summary()
print(result.credible_effects())
```

`set_fdr(est_fdr=0.1)` chooses the spike-and-slab threshold for the desired expected FDR; credible effects are the populations whose change is supported relative to the reference.

## sccomp - outlier-robust Bayesian composition (R)

**Goal:** Test composition (and variability) jointly, robust to outlier samples.

**Approach:** Estimate the beta-binomial model from a count table or cell-level data with `sccomp_estimate`, optionally remove outliers, then test contrasts with `sccomp_test`, which returns a Bayesian FDR (`c_FDR`).

```r
library(sccomp)

res <- counts_tbl |>
    sccomp_estimate(formula_composition = ~ condition, .sample = sample, .cell_group = cell_type, .count = count, cores = 1) |>
    sccomp_remove_outliers(cores = 1) |>
    sccomp_test()
res[res$c_FDR < 0.05, c('cell_type', 'c_effect', 'c_FDR')]
```

The block uses the sccomp 1.x tidy-eval arguments (`.sample`, `.cell_group`, `.count`, bare column names). The current upstream README (checked, not run) writes them as strings without the dot (`sample = "sample", cell_group = "cell_type", abundance = "count"`) and names the result column `cell_group`; check `?sccomp_estimate` on the installed version and adapt. sccomp needs CmdStan, which was not installed in the audit env, so this block has never been executed here.

`sccomp_test` reports `c_effect` (composition log-fold change) and `c_FDR`; modeling variability separately catches groups that differ in dispersion, not just mean proportion.

## propeller - fast frequentist proportions (R)

**Goal:** Quickly test cell-type proportion differences across groups.

**Approach:** Compute per-sample proportions, apply a logit (speckle 1.6.0 default; `transform = 'asin'` gives arcsin-sqrt) variance-stabilizing transform, and run a limma moderated test per cell type.

```r
library(speckle)

out <- propeller(clusters = seurat_obj$cell_type, sample = seurat_obj$sample, group = seurat_obj$condition)
out[out$FDR < 0.05, ]
```

Audit run (speckle 1.6.0): the default call prints "Performing logit transformation of proportions" and called only NK on the 8-donor set (planted truth: NK only). `transform = 'asin'` is more liberal on the same data: NK FDR 1e-9 but also CD8 T cells at FDR 0.009, a false positive, so keep the logit default unless a method needs asin. With one sample per group it aborts with `.ebayes(...): No finite residual standard deviations` (see Common Errors). Recovery by group size on the same data (NK, planted expansion): n=2/group FDR 0.012, n=3 FDR 6e-4, n=4 FDR 9e-5, and no other cell type was called at any n.

propeller is the fast default for several samples per group; for outliers, over-dispersion, or random effects, prefer sccomp or scCODA.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Many cell types flagged as changed, all anti-correlated | Per-cluster proportion t-tests ignore the simplex | Use scCODA/sccomp/propeller/Milo, which model the joint composition |
| scCODA verdict flips between runs | Reference cell type mis-chosen or actually changing | Pick a stable abundant reference, or `reference_cell_type='automatic'` |
| No significant abundance change despite an obvious shift | Too few biological replicates; underpowered | Add donors (not cells); report effect sizes / credible intervals |
| Milo returns zero significant neighborhoods (raw p also > 0.05) while a cluster-level test is highly significant | Neighborhood GLM under-powered at this cell count (see the Milo power note) | Report the cluster-level result as primary and Milo's logFC as effect size; do not read the Milo null as no change |
| propeller aborts with `No finite residual standard deviations` | One sample per group: no residual degrees of freedom | The effect is unidentifiable; add biological replicates (>=2, ideally 3-4 per group), do not switch tools |
| Different neighborhoods or a different posterior on each run | `makeNhoods` sampling and scCODA HMC are unseeded | Set `set.seed()` / `tf.random.set_seed()` and report seed plus k/`prop` or HMC settings |
| Milo neighborhoods look noisy / unstable | k or `prop` too small, or embedding not integrated | Increase k/prop; build the graph on a batch-corrected reduced dim |
| "DE genes" between conditions but expression unchanged per cell | Compositional shift masquerading as DE | Run a differential-abundance test alongside the DE analysis |
| propeller p-values too liberal with few samples | Frequentist test under-powered/over-confident at small n | Use a Bayesian model (sccomp/scCODA) and report uncertainty |
| Abundance significant only in one direction across all types | Reporting raw proportions without the constraint | Interpret relative to a reference and report which population actually drives the shift |

## Related Skills

- clustering - Define the clusters whose abundance is tested (cluster-based methods)
- cell-annotation - Annotate cell types before testing their proportions
- markers-annotation - Pair condition DE with abundance testing to separate the confound
- batch-integration - Build the integrated embedding Milo's kNN graph relies on
- differential-expression/deseq2-basics - Pseudobulk condition DE that abundance testing complements
- pathway-analysis/go-enrichment - Characterize the populations that expanded or contracted

## References

- Dann et al. 2022, Nat Biotechnol 40:245-253 - Milo; differential abundance on kNN-graph neighborhoods with SpatialFDR.
- Buttner et al. 2021, Nat Commun 12:6876 - scCODA; Bayesian Dirichlet-multinomial compositional analysis with a reference cell type.
- Mangiola et al. 2023, PNAS 120(33):e2203828120 - sccomp; outlier-robust Bayesian differential composition and variability.
- Phipson et al. 2022, Bioinformatics 38(20):4720 - propeller; arcsin-sqrt transform plus limma for cell-type proportion testing.
- Squair et al. 2021, Nat Commun 12:5692 - sample, not cell, is the unit of replication for cross-condition single-cell claims.
