---
name: bio-workflows-proteomics-pipeline
category: Data Analysis
description: Orchestrates bottom-up proteomics from a search engine's output (MaxQuant/DIA-NN) to differential protein abundance with limma/DEqMS/MSstats. Use when committing the search database + acquisition mode (DDA vs DIA) up front, re-controlling FDR at PSM AND peptide AND protein-group level (not just PSM), removing contaminant/reverse rows and inspecting RAW distributions before normalizing, bridging cross-plex TMT with an IRS reference channel, modeling MNAR missingness rather than downshift-imputing on/off proteins, batching as a covariate (not pre-subtracted), and testing with treat()/DEqMS. Hands mechanism to the proteomics component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: limma
workflow: true
depends_on:
  - proteomics/data-import
  - proteomics/proteomics-qc
  - proteomics/quantification
  - proteomics/protein-inference
  - proteomics/differential-abundance
  - proteomics/dia-analysis
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: MSnbase 2.32.0, limma 3.62.2, MSstats 4.14.2, DEqMS 1.24.0, proDA 1.20.0, MSstatsTMT 2.14.2, arrow 23.0.1 (DIA-NN report.parquet), dplyr 1.2.1, tidyr 1.3.2, ggplot2 4.0.3 (checked 2026-09-15 on R 4.4.3)

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Proteomics Pipeline

**"Process my proteomics data from raw MS files to differential abundance"** -> Orchestrate data import (pyopenms/MaxQuant), QC assessment, protein quantification, normalization, differential abundance testing (limma/DEqMS, or MSstats for feature-level designs), and PTM analysis.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

Scope: research cohort comparison. A pipeline output is not a diagnostic -- do not use it to triage, escalate or choose treatment for an individual patient; route those questions to validated clinical assays and the treating clinician.

## The governing principle

Bottom-up proteomics never measures proteins; it measures peptides and INFERS proteins, and the trustworthiness decisions are made at seams before the statistics.

1. **The search database + acquisition mode are committed once and inherited by everything.** The FASTA fixes the target-decoy frame (concatenated one-search FDR = #decoy/#target; a separate-search design needs mix-max instead — mixing the two mis-estimates FDR), what counts as a "unique peptide" (relative to the DB: canonical vs +isoforms), and the contaminants (cRAP must be IN the search DB from the start; a contaminant can BE the protein of interest, so never blind-delete `CON__` rows). DDA vs DIA is set at the instrument and dictates which imputation is even legitimate.
2. **FDR is re-controlled at THREE levels, not just PSM — and match-between-runs has its OWN FDR.** 1% PSM-FDR does NOT give 1% protein-FDR — each level (PSM, peptide, protein-group) needs its own target-decoy estimation; PSM-only filtering yields 10-30% real protein-FDR on deep data (one false PSM nucleates a false one-hit-wonder, and false proteins grow with dataset size). Use picked-protein/picked-group FDR. The two-peptide rule INCREASES protein-FDR, it does not reduce it. MBR transfers IDs across runs by RT/m-z and can be wrong for low-abundance precursors — do NOT report MBR-filled counts as directly measured; DIA-NN controls MBR-FDR via `Lib.*` q-values, IonQuant via an explicit MBR-FDR mixture model.
3. **Missingness is MODELED, not filled.** DDA missingness is structured left-censored MNAR; downshift imputation (mean=mu-1.8sigma) on an on/off protein inflates the t-numerator AND deflates the denominator (the volcano "wing" artifact). The honest report for a protein missing in one whole group is "undetected in group B", not a fold change — model the MNAR (proDA/msqrob2/MSstats-AFT).
4. **Normalize AFTER contaminant removal and AFTER inspecting raw distributions; batch is a covariate, not pre-subtracted.** Median-normalizing first mathematically erases a 3x-low load. Cross-plex TMT is invalid without an IRS bridge. `removeBatchEffect` before testing understates residual variance (anticonservative p) — put batch in the same model.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Search FASTA + target-decoy strategy | The FDR estimator, what a "unique peptide" is, which contaminants exist; a mismatch mis-estimates FDR silently |
| Enzyme + fixed/variable mods (Carbamidomethyl-Cys fixed) | Which peptides exist to quantify; a fixed-mod misconfig loses all Cys peptides |
| DDA vs DIA acquisition mode | Missingness structure (MNAR vs ~MCAR), whether TMT is possible, which imputation is legitimate |
| FDR framing (PSM + peptide + protein-group, 1% each) | Real protein-FDR; PSM-only is 10-30% wrong on deep data |

## Pipeline Overview

```
Raw MS Data (mzML) --> MaxQuant/DIA-NN --> proteinGroups.txt
                                                 |
                                                 v
            +--------------------------------------------+
            |             proteomics-pipeline            |
            +--------------------------------------------+
            |  1. Data Import & Filtering                |
            |  2. Log2 + inspect RAW distributions       |
            |  3. Normalization (after the inspection)   |
            |  4. Per-Group Completeness Filter          |
            |  5. QC: PCA, Correlation                   |
            |  6. Differential Abundance (limma/MSstats) |
            |  7. Visualization & Export                 |
            +--------------------------------------------+
                                                 |
                                                 v
                  Differential Proteins + Volcano Plots
```

## Inputs and Install

- MaxQuant: `proteinGroups.txt` for the limma route; `evidence.txt` + `proteinGroups.txt` + `annotation.csv` (MSstats route). DIA-NN: `report.parquet`.
- `sample_annotation.csv` (limma route), one row per sample: `sample` (must equal the intensity column name), `condition`, `replicate`, `batch`. `batch` is required whenever samples were acquired in more than one run/day/plex -- the design branch keys on that column and puts batch in the model as a covariate; without it the batch effect stays in the residual. `condition` may have more than two levels (dose series, time course); every non-reference level is contrasted against the first level.

```csv
sample,condition,replicate,batch
Sample1,Control,1,B1
Sample2,Control,2,B2
Sample3,Treatment,1,B1
Sample4,Treatment,2,B2
```

- Install: `BiocManager::install(c('limma', 'DEqMS', 'proDA', 'MSstats', 'MSstatsTMT', 'MSnbase'))`; `install.packages(c('pheatmap', 'ggplot2', 'arrow', 'dplyr', 'tidyr'))`.

## Complete R Workflow

**Goal:** Turn a MaxQuant or DIA-NN protein matrix into a table of differentially abundant proteins with honest missing-value handling.

**Approach:** Strip bookkeeping rows, log2 and inspect the RAW per-sample distributions (dropping failed loads before normalization can hide them), median-center, filter on per-group completeness, then test the OBSERVED values with moderated limma using treat() for a minimum fold change and batch as a covariate -- nothing is imputed. Upgrade to proDA when the dropout itself has to be modeled.

TMT/iTRAQ, SILAC, DIA-NN `report.parquet` and MSstats feature-level input have their own routes: see Reference Files below.

```bash
Rscript scripts/limma_pipeline.R proteinGroups.txt sample_annotation.csv proteomics_results.csv [min_frac=0.6] [min_fold_change=1.5]
```

`scripts/limma_pipeline.R` runs the steps below in order; its comments carry the reasoning, so read them before changing a threshold.

1. Import with `quote = ''` and `comment.char = ''`, check the row count against `readLines`, drop contaminant / reverse / only-by-site rows with `%in% '+'`.
2. log2 (0 -> NA), then inspect the RAW `Intensity.` columns (not LFQ): drop a sample whose raw median is >= 1 log2 low, or whose ID count is more than 3 MADs below the median.
3. Median-normalize (assumes most proteins are unchanged and the changes roughly symmetric; see Common Errors).
4. Per-group completeness filter (>= `min_frac` present in at least one group). Nothing is imputed; the proDA upgrade and the downshift DO NOT are in the script's comments.
5. QC: PCA on complete cases, guarded, with a pairwise-correlation fallback.
6. limma on the observed values: `~ 0 + condition` plus `factor(batch)` when batch has more than one level, contrasts built from the condition levels against the first, `treat(lfc = log2(min_fold_change), trend = TRUE, robust = TRUE)`, one global BH over every protein x contrast. `min_fold_change` is a floor on EVERY contrast: on a dose series or time course it can zero the intermediate levels (measured, 4 replicates per level, 120 proteins with a real 0.7 log2 low-dose effect: High / Low calls 20 / 0 at log2(1.5), 52 / 0 at lfc 0.3, 78 / 8 at lfc 0), so lower it and screen with the eBayes F-test given in the script's comments.
7. Output: one row per protein x contrast (`protein`, `contrast`, `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `significant`). `significant` is the call to use (global decideTests, not `adj.P.Val < 0.05` alone); NA `logFC` = not estimable, i.e. undetected in one group.

## QC Checkpoints

| Stage | Check | Action if Failed |
|-------|-------|------------------|
| Import | >1000 proteins | Re-run MaxQuant |
| Filter | <30% removed | Check sample prep |
| Missing | <40% per sample | Check MS performance |
| PCA | Replicates cluster | Check for batch effects |
| Design | >= 3 biological replicates per condition | Do not run a per-protein test; report as exploratory or add replicates |
| Stats | FC/FDR pre-specified | Verify thresholds were pre-specified; inspect the volcano for downshift-imputation 'anchor arms' |

## Reference Files

The main path above (limma on a MaxQuant `proteinGroups.txt` matrix, `scripts/limma_pipeline.R`) is the default. Read a reference file and run its script only when the request needs its route:

| Route | Read | Run | When |
|-------|------|-----|------|
| MSstats feature-level model | `references/msstats.md` | `scripts/msstats_maxquant.R` | peptide/feature-level input (`evidence.txt`), feature-level mixed models, contrast matrices for 3+ conditions |
| TMT / iTRAQ | `references/tmt-isobaric.md` | `scripts/tmt_impurity_correct.R` (single plex), `scripts/msstatstmt_multiplex.R` (2+ plexes) | isobaric reporter data; CoA impurity correction; multi-plex experiments (reference-channel bridge via MSstatsTMT) |
| SILAC | `references/silac.md` | `scripts/silac_limma.R` | MaxQuant SILAC H/L ratios |
| DIA-NN | `references/dia-nn.md` | `scripts/diann_matrix.R` | DIA-NN `report.parquet` instead of `proteinGroups.txt` |

Every script takes its inputs as arguments (usage line in its header).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| ~1% PSM-FDR but 10-30% wrong proteins | FDR controlled only at PSM level | Estimate FDR at peptide AND protein-group level (picked-group FDR) |
| Volcano "wings" of huge-FC on/off proteins | Downshift imputation on MNAR (Perseus/MaxQuant) | Model the MNAR (proDA/msqrob2/MSstats-AFT); report "undetected in group B", not a fold change |
| Cross-plex TMT ratios differ 2-5x for no biology | Compared TMT across plexes without IRS | Pooled reference channel in EVERY plex + IRS bridge before comparison |
| A failed-load sample silently carried forward | Normalized before inspecting raw distributions | Filter contaminant/reverse rows -> inspect raw boxplots + ID counts -> remove failures -> THEN normalize |
| Anticonservative p-values | `removeBatchEffect` before testing | Put batch in the model (`~ batch + condition`); removeBatchEffect only for PCA |
| Every ratio subtly wrong | Wrong intensity column (`Intensity` vs `LFQ intensity` vs `iBAQ`) | Pick the right column; convert 0 -> NaN before log2 |
| Spurious DA that flips between conditions | Razor-peptide inference reassigns a shared peptide | Quantify at protein-group level or unique-peptides-only for sensitive comparisons |
| `EOF within quoted string`; far fewer rows than the file has lines | default `read.table` quoting on MaxQuant tables with apostrophes (`5'-nucleotidase`) and `#` | `read.table(..., quote = '', comment.char = '')`, then check `nrow` against `readLines` |
| `eBayes`/`lmFit`: `missing value where TRUE/FALSE needed` on a DIA-NN matrix | DIA-NN writes 0 for "not quantified"; `log2(0)` is `-Inf` | `m[m == 0] <- NA` before `log2`, then `stopifnot(!any(is.infinite(m)))` |
| `t.test`: `not enough 'x' observations` partway through a SILAC run | a protein quantified in one replicate only | filter to >= 2 finite ratios; moderated one-sample limma; report BH-adjusted p, never raw |
| Impurity correction makes adjacent TMT10 channels worse, not better | `makeImpuritiesMatrix(filename=)` places CoA Da-offsets by POSITION, but TMT10/TMTpro interleave N and C | build the matrix by channel name (rows = source reagent) and pass it straight to `purityCorrect` |
| A 3x-low injection sails through the raw-distribution check | the check read `LFQ intensity`, which MaxLFQ already renormalized, and used a 50%-of-median ID rule | inspect the raw `Intensity.` columns; flag a >= 1 log2 load shift or an ID count > 3 MADs below the median |
| `prcomp`: `a dimension is zero` at the QC step, before any statistics | no protein is observed in EVERY sample, so the complete-case matrix is empty (12 samples at realistic MNAR dropout is enough) | guard `nrow(complete) >= 3`; fall back to `cor(..., use = 'pairwise.complete.obs')`, lower `min_frac`, or drop the sparsest samples -- never impute to fill the PCA matrix |
| `makeContrasts`: `object 'Treatment' not found`, with a correct design matrix one line above | the contrast is hard-coded for two conditions; a dose series or time course has three or more | build the contrasts from `levels(sample_info$condition)` against the reference level, and adjust ACROSS them with `decideTests(method = 'global')` |
| Every unchanged protein drifts one way; the hit list is implausibly one-directional | median centering / `equalizeMedians` on a design whose changes are NOT symmetric (pulldown, secretome, strong one-sided response) | normalize on a set expected to be unchanged (spike-ins, `globalStandards`); check that the mean log2FC over expected-null proteins is ~0 |
| Impurity correction still leaves adjacent-channel bleed, no error, no negative values | the lot CoA was transposed -- same names, same shape, zero negatives, 2.7x worse than the right orientation | check orientation by row vs column sums to 100% before `purityCorrect`; rows are the SOURCE reagent |

## References

- Elias JE, Gygi SP (2007) Target-decoy search strategy for increased confidence in large-scale protein identifications by mass spectrometry. *Nature Methods* 4:207-214. DOI 10.1038/nmeth1019.
- Savitski MM, Wilhelm M, Hahne H, Kuster B, Bantscheff M (2015) A scalable approach for protein false discovery rate estimation in large proteomic data sets. *Molecular & Cellular Proteomics* 14:2394-2404. DOI 10.1074/mcp.M114.046995. (picked-protein FDR.)
- Plubell DL, Wilmarth PA, Zhao Y, et al (2017) Extended multiplexing of tandem mass tags (TMT) labeling reveals age and high-fat-diet specific proteome changes in mouse epididymal adipose tissue. *Molecular & Cellular Proteomics* 16:873-890. DOI 10.1074/mcp.M116.065524. (IRS.)
- Ritchie ME, Phipson B, Wu D, et al (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 43:e47. DOI 10.1093/nar/gkv007.
- Zhu Y, Orre LM, Zhou Tran Y, et al (2020) DEqMS: a method for accurate variance estimation in differential protein expression analysis. *Molecular & Cellular Proteomics* 19:1047-1057. DOI 10.1074/mcp.TIR119.001646.

## Related Skills

- proteomics/data-import - Load MS data formats
- proteomics/proteomics-qc - Quality control before analysis
- proteomics/quantification - Normalization, TMT IRS bridge, SILAC mechanics
- proteomics/protein-inference - Razor/shared-peptide assignment to protein groups
- proteomics/differential-abundance - Modeling missingness, moderated testing details
- proteomics/dia-analysis - DIA-NN report parsing and q-value filtering
- proteomics/ptm-analysis - Phosphoproteomics and other PTMs
- data-visualization/volcano-and-ma-plots - Volcano plots with LFC shrinkage
