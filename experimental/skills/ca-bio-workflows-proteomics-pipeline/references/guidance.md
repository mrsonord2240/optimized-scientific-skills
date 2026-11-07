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


## QC Checkpoints

| Stage | Check | Action if Failed |
|-------|-------|------------------|
| Import | >1000 proteins | Re-run MaxQuant |
| Filter | <30% removed | Check sample prep |
| Missing | <40% per sample | Check MS performance |
| PCA | Replicates cluster | Check for batch effects |
| Design | >= 3 biological replicates per condition | Do not run a per-protein test; report as exploratory or add replicates |
| Stats | FC/FDR pre-specified | Verify thresholds were pre-specified; inspect the volcano for downshift-imputation 'anchor arms' |


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
