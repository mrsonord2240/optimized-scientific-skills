---
name: bio-workflows-metabolomics-pipeline
description: Orchestrates the untargeted LC-MS metabolomics pipeline end-to-end (xcms 4.x feature extraction, QC/drift/normalization, confidence-stratified annotation, permutation-validated statistics, background-aware pathway mapping), naming what each stage decides and where it silently fails. Use when running a full LC-MS metabolomics study from raw mzML to enriched pathways and needing the honest handoffs between stages. Each stage defers to its component skill for parameters and traps; for stable-isotope flux (a separate branch, not this untargeted flow) see metabolomics/isotope-tracing.
tool_type: r
primary_tool: xcms
workflow: true
depends_on:
  - metabolomics/xcms-preprocessing
  - metabolomics/metabolite-annotation
  - metabolomics/normalization-qc
  - metabolomics/statistical-analysis
  - metabolomics/pathway-mapping
  - metabolomics/lipidomics
  - metabolomics/targeted-analysis
  - metabolomics/msdial-preprocessing
qc_checkpoints:
  - after_extraction: "Feature count plausible after redundancy collapse; EICs of top hits inspect cleanly; is_filled cells tracked"
  - after_drift: "QC RSD DROPS after correction AND biological-sample RSD is unchanged (a rise means the spline absorbed signal)"
  - after_qc_filter: "QC RSD <=20-30%, D-ratio <=0.5, blank ratio >=3-5x, detection rate >=50-80% applied BEFORE imputation"
  - before_stats: "No NAs remain in the matrix handed to OPLS-DA/PLS-DA (impute by mechanism first; opls() errors on NAs, does not warn)"
  - after_stats: "Univariate FDR (BH) AND permutation-validated multivariate (permI>=1000; Q2 high with small pQ2); getSummaryDF() has >0 rows, not a silent empty model"
  - after_annotation: "MSI/Schymanski level assigned per compound; only Level 1-2 enter identified-ORA"
  - after_pathway: "Background = assay coverage (identified ORA) OR full feature table (mummichog); PREDICTED vs MEASURED stated"
license: MIT
---

## Version Compatibility

Reference examples tested with: xcms 4.x+ (MsExperiment/XcmsExperiment), pmp 1.14+, ropls 1.34+, MetaboAnalystR 4.0+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

Checked on pmp 1.18.0, ropls 1.38.0, MetaboAnalystR 4.3.0. Stage 4's `fit_discriminant_guarded()`
and Stage 5's `current.msg`/`Setup.KEGGReferenceMetabolome`/Local-Only-ORA pattern are defined in
metabolomics/statistical-analysis and metabolomics/pathway-mapping respectively -- this Skill
calls them by name rather than redefining them, so a fix to either component Skill's function
applies here too.

This pipeline is only as honest as its weakest stage: a flawless feature table fed to a too-flexible drift model, or a clean OPLS-DA plot fed to background-free enrichment, produces confident wrong biology. Validate each stage against its own held-out check (QCs, permutation null, assay-coverage background), not against the next stage looking nice.

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metabolomics Pipeline

**"Process my LC-MS metabolomics data end-to-end"** -> Chain xcms feature extraction, QC/normalization, confidence-stratified annotation, validated statistics, and background-aware pathway mapping, treating each stage's output as a hypothesis its component skill scrutinizes.
- R: `readMsExperiment()` -> `findChromPeaks()` -> `groupChromPeaks()` -> `fillChromPeaks()` -> `featureValues()` (xcms), then `QCRSC()`/`pqn_normalisation()` (pmp), `opls()` (ropls), `CalculateOraScore()`/`PerformPSEA()` (MetaboAnalystR)

## The governing principle

An untargeted metabolomics result is only as honest as its weakest seam: a feature table is a PARAMETERIZED HYPOTHESIS, not a measurement; injection-order drift is technical variance collinear with run order that masquerades as biology unless corrected BEFORE stats; and a database name is an MSI Level 4-5 guess until an authentic standard makes it Level 1. Three commitments are fixed at the bench and inherited by everything downstream:

1. **Ionization mode is a mode-lock.** Positive and negative ESI have entirely different adduct chemistries ([M+H]+/[M+Na]+ vs [M-H]-/[M+Cl]-); the mode selects which mass-shift table is legal for every candidate, so annotation and pathway mapping must run against the SAME mode the feature was acquired in, and mixed-mode data must carry a per-feature mode column all the way through.
2. **The annotation-confidence contract gates entry into pathway mapping.** MSI/Schymanski level is a made-once reporting decision: a bare DB hit with no orthogonal MS/MS or RT evidence is Level 4-5, NOT Level 2-3, and feeding tentative IDs into ORA as if confirmed launders uncertainty into a pathway p-value.
3. **The pooled-QC + blank + dilution baseline, fixed at the bench, makes every correction possible.** The pooled QC is the substrate for BOTH drift correction and feature-quality filtering; biological samples MUST be block-randomized in run order (group confounded with order is the unwinnable "original sin"); conditioning injections are excluded from drift modeling.

## What Each Stage Decides and Where the Traps Are

This skill is an orchestrator: it sequences the five component skills and enforces the honest handoffs between them. It does not re-teach each stage's parameters -- those live in the component SKILLs cited per row.

| Stage | The decision it owns | The trap it must not paper over | Defers to |
|---|---|---|---|
| 1. Feature extraction | centWave/grouping/alignment parameters that set the detection floor | A feature table is a parameterized hypothesis; `fillChromPeaks` fabricates intensities; 1 compound = 5-15 features | metabolomics/xcms-preprocessing |
| 2. QC + normalization | drift correction, RSD/D-ratio filtering, dilution normalization, mechanism-aware imputation | Over-correction is invisible to QC RSD; half-min-impute-then-test inflates significance; confounded design is unrescuable | metabolomics/normalization-qc |
| 3. Annotation | the MSI/Schymanski confidence level of every name | A database hit is Level 4-5, not an identification; ambiguous m/z inflates downstream pathways | metabolomics/metabolite-annotation |
| 4. Statistics | univariate FDR + permutation-validated multivariate, reconciled | A clean PLS-DA score plot is the generic output of p>>n; R2 is no evidence; scaling changes conclusions | metabolomics/statistical-analysis |
| 5. Pathway mapping | ORA on IDs vs mummichog on m/z, with an explicit background | The background IS the null; enrichment launders annotation uncertainty into confident biology | metabolomics/pathway-mapping |

## Pipeline Flow

```
raw mzML (centroided)
   |  metabolomics/xcms-preprocessing
   v  readMsExperiment -> findChromPeaks -> adjustRtime -> groupChromPeaks -> fillChromPeaks
features x samples table (+ is_filled flags, mzmed/rtmed)
   |  metabolomics/normalization-qc
   v  blank/detection filter -> within-batch drift (QCRSC) -> RSD/D-ratio filter -> PQN -> mechanism-aware impute
QC-clean, dilution-normalized matrix
   |  split: statistics  AND  annotation (independent axes)
   v
metabolomics/statistical-analysis            metabolomics/metabolite-annotation
permutation-validated hits + univariate FDR  confidence-stratified names (MSI level per feature)
   |                                                |
   +-------------------- join on feature_id --------+
   v  metabolomics/pathway-mapping
identified compounds -> ORA/MSEA   OR   raw m/z (no IDs) -> mummichog/PSEA (background = FULL table)
   v
pathways "consistent with perturbation", conditional on annotations + background
```

Stable-isotope tracing (flux) is a SEPARATE branch off labeled raw data, not a stage of this untargeted flow -- see metabolomics/isotope-tracing.

## Stage 1 -- Feature Extraction (modern xcms 4.x)

**Goal:** Turn centroided mzML into a features-by-samples table, carrying the parameters as part of the result.

**Approach:** Use the `MsExperiment`/`XcmsExperiment` containers with `*Param` objects; align to pooled QC, group AFTER alignment (obiwarp aligns the raw profile directly, so no pre-grouping is needed; the PeakGroups method instead needs group -> align -> regroup because it uses grouped anchor peaks), and treat filled values as imputations. Full parameter rationale (ppm, peakwidth, bw, prefilter) lives in metabolomics/xcms-preprocessing.

```r
library(xcms)
# pd: data.frame, one row per file, with a sample_group column ('QC'/'Control'/'Treatment')
# ionization_mode: 'positive' or 'negative', set from the acquisition method -- carried through
# to Stage 3/5 as defs$mode below (commitment #1, the mode-lock). A mixed-mode study runs this
# whole stage twice, once per mode, and merges the resulting feature tables afterward.
raw <- readMsExperiment(spectraFiles = mzml_files, sampleData = pd)

cwp <- CentWaveParam(ppm = 10, peakwidth = c(2, 20), snthresh = 10,
                     prefilter = c(3, 1000), noise = 1000)   # set from instrument; see xcms-preprocessing
xdata <- findChromPeaks(raw, param = cwp)
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6,
    subset = which(sampleData(xdata)$sample_group == 'QC'), subsetAdjust = 'average'))   # anchor RT alignment on pooled QCs
pdp <- PeakDensityParam(sampleGroups = sampleData(xdata)$sample_group,
                        bw = 5, minFraction = 0.5, binSize = 0.025)
xdata <- groupChromPeaks(xdata, param = pdp)        # group on corrected RT (obiwarp needs no pre-grouping)
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())

feat <- featureValues(xdata, value = 'into')        # features x samples; filled cells are imputations
defs <- featureDefinitions(xdata)                   # mzmed / rtmed per feature, for annotation + mummichog
```

## Stage 2 -- QC, Drift, Normalization (not naive median + half-min)

**Goal:** Filter junk features, flatten injection-order drift, normalize per-sample dilution, and impute by mechanism -- before any test sees the data.

**Approach:** Follow the normalization-qc pipeline order: blank/detection filter -> within-batch drift correction (QCRSC) -> RSD/D-ratio filter -> PQN -> mechanism-aware imputation. Do NOT silently half-min-impute and feed limma; validate drift correction on held-out QCs, not on QC clustering.

```r
library(pmp)
library(imputeLCMD)
# feat (Stage 1 featureValues() output) is already features x samples -- the pmp convention.
# No transpose needed; assert it instead of assuming, so an API change is caught, not silently
# masked (pmp's own check_peak_matrix re-transposes with no warning when it can, which is what
# hid this exact bug before: a wrong `t(feat)` here still "worked").
fm <- feat
stopifnot(nrow(fm) == nrow(defs), ncol(fm) == length(sample_class))

filtered <- filter_peaks_by_fraction(fm, classes = sample_class, min_frac = 0.5, qc_label = 'QC')
corrected <- QCRSC(df = filtered, order = injection_order, batch = batch_id,
                   classes = sample_class, spar = 0, minQC = 5, qc_label = 'QC')  # CV-selected spline
rsd_filtered <- filter_peaks_by_rsd(corrected, max_rsd = 30, classes = sample_class, qc_label = 'QC')
normalized <- pqn_normalisation(rsd_filtered, classes = sample_class, qc_label = 'QC')

# QCRSC silently returns an ENTIRE batch as all-NA when that batch has fewer QCs than minQC (see
# normalization-qc's Common Errors) -- verified on real MTBLS79 data: 5 of 8 batches had only 4
# QCs (< minQC=5), wiping 90 of 172 samples. A wholly-missing sample carries NO measured
# information to impute -- report and drop it; do not fabricate a profile for it.
nm <- as.matrix(normalized)
wiped <- colMeans(is.na(nm)) == 1
if (any(wiped)) {
  cat(sum(wiped), 'of', ncol(nm), 'samples came back all-NA (QCRSC: their batch had < minQC QCs)',
      '-- dropping, not imputing:\n')
  print(table(batch_id[wiped]))
  nm <- nm[, !wiped, drop = FALSE]
  sample_class <- sample_class[!wiped]   # keep every per-sample vector in sync for Stage 4
}

# The remaining holes are the sparse, mechanism-driven kind imputation is actually valid for.
# QRILC (MNAR / left-censored) is only valid on log-scale intensities; on raw intensities it
# silently draws negative (impossible) values. Round-trip through log2/2^x, per normalization-qc.
log_mat <- log2(nm)
imputed <- 2^(impute.QRILC(log_mat, tune.sigma = 1)[[1]])
stopifnot(min(imputed, na.rm = TRUE) >= 0, !anyNA(imputed))  # no NAs may reach Stage 4 (opls() cannot tolerate them)
```

Drift correction should lower QC RSD AND leave biological-sample RSD unchanged; if biological RSD rises, the spline absorbed signal. Mechanism-aware imputation (QRILC/GSimp for left-censored zeros) replaces the old half-min step, which collapses imputed-subset variance and inflates false significance. `imputed` (not `normalized`) is what Stage 4 receives -- feeding `normalized` directly, unimputed and still carrying any wiped samples, is exactly the bug that crashes real multi-batch data (see Stage 4 below). Verified end to end on real MTBLS79 data (2433 features, 172 samples, 8 batches): 90 wiped samples dropped, 82 survive with only sparse residual NAs (max 18% per sample), QRILC imputes cleanly, and the Stage 4 OPLS-DA below fits successfully (`pR2Y = pQ2 = 0.001`).

## Stage 3 -- Annotation Before Claiming IDs

**Goal:** Attach an MSI/Schymanski confidence level to each feature so the pathway stage knows what it is allowed to claim.

**Approach:** Match MS/MS to a library (Level 2a) or run SIRIUS/CSI:FingerID (formula Level 4, structure Level 2b/3); a bare m/z is Level 5. Collapse ion families (CAMERA) first so adducts of one compound are not counted as separate metabolites. Mechanics and thresholds live in metabolomics/metabolite-annotation. Annotation and statistics are independent axes -- run them in parallel and join on feature_id.

```r
# Mode-lock (commitment #1): every feature carries the ionization mode it was acquired in;
# annotation and pathway mapping must use the SAME mode's adduct/mass-shift table for that
# feature. Thread the column through explicitly -- do not assume single-mode just because
# most studies are. 'positive'/'negative' is set per acquisition, never inferred from m/z.
defs$mode <- ionization_mode
stopifnot(!anyNA(defs$mode))
if (length(unique(defs$mode)) > 1) {
  # Mixed-mode study (separate pos/neg acquisitions merged before this stage): split, annotate
  # each half against ITS OWN adduct table (metabolite-annotation), then recombine. Never score
  # a positive-mode candidate against a negative-mode feature or vice versa.
  by_mode <- split(defs, defs$mode)
}
```

## Stage 4 -- Statistics (univariate FDR + validated multivariate)

**Goal:** Decide which metabolites genuinely differ, with neither a score plot nor an unadjusted p-value standing alone.

**Approach:** Transform (if heteroscedastic), pick a scaling explicitly (run >=1 alternative and check the conclusion is not scaling-fragile), run a Welch/Mann-Whitney univariate test with BH FDR, AND a permutation-validated OPLS-DA (`permI >= 1000`), then reconcile the two. Full validation checklist in metabolomics/statistical-analysis.

```r
library(ropls)
# t(imputed): samples x features; sample_class carried over from Stage 2 (post drop/impute) so
# group labels stay aligned to imputed's columns. Use `imputed`, NOT `normalized` -- opls()
# cannot tolerate NAs, and a real multi-batch run is not sparse enough to skip this.
study_samples <- sample_class != 'QC'
group <- factor(sample_class[study_samples])
x <- t(imputed)[study_samples, ]
stopifnot(!anyNA(x))   # last check before the hand-off; see Stage 2's imputation step above

# ropls's own CV significance test on the first predictive component can reject it and silently
# return a 0-row summaryDF / empty model with no warning under info.txtC='none' (~40% of runs in
# a p>>n regime) -- use the guarded fit from metabolomics/statistical-analysis, not a bare opls()
# call, so this hand-off cannot silently produce an unusable model:
fit <- fit_discriminant_guarded(x, group, scaleC = 'pareto', permI = 1000)   # see statistical-analysis
oplsda <- fit$model
cat('Model type actually fit:', fit$type, '\n')   # report this -- it is not always OPLS-DA
summ <- getSummaryDF(oplsda)   # claim licensed only if Q2 high AND pQ2 small; R2Y alone proves nothing
```

Univariate Welch + BH (`p.adjust(method='BH')` in R, `multipletests(method='fdr_bh')` in Python -- neither default is BH) gives the per-feature answer with effect sizes. Features are correlated (adducts, pathways), so collapse to compounds before counting "how many metabolites changed."

## Stage 5 -- Pathway Mapping (the background is the null)

**Goal:** Interpret the differential result in pathway context without laundering annotation uncertainty into confident biology.

**Approach:** Two disjoint entry points. Confidently identified compounds -> ORA/MSEA with an assay-coverage background (NOT all of KEGG). Raw m/z with no IDs -> mummichog/PSEA whose permutation null is sampled from the FULL feature table. Either way, report mapping coverage and the MSI levels of the driving compounds; downgrade claims to "consistent with perturbation." Full method choice and background construction in metabolomics/pathway-mapping.

```r
# MSI/Schymanski gate (commitment #2): only Level 1-2 (1, 2a, 2b) may enter identified-ORA as an
# "identification" -- Level 3-5 is a database-name hypothesis, not a name (see
# metabolite-annotation). This filter is the one place the principle is enforced as code, not
# just stated in a QC-checkpoint row.
# annotated: Stage 3's per-feature output (feature_id, name, msi_level -- msi_level values like
# 1, '2a', '2b', 3, 4, 5 per metabolite-annotation's assign_level()), joined to feature_id.
identified_compounds <- unique(annotated$name[grepl('^[12]', as.character(annotated$msi_level))])
stopifnot(length(identified_compounds) > 0)   # all Level 3-5? use Path B (mummichog) instead, not a forced ORA

current.msg <- character(0); err.vec <- character(0)  # required outside the Shiny app; see pathway-mapping
library(MetaboAnalystR)
# Path A: identified compounds (MSI level 1-2) -> ORA
mSet <- InitDataObjects('conc', 'pathora', FALSE)
mSet <- SetOrganism(mSet, 'hsa')
mSet <- Setup.MapData(mSet, identified_compounds)
mSet <- CrossReferencing(mSet, 'name')
mSet <- CreateMappingResultTable(mSet)              # inspect coverage before trusting any p-value
mSet <- SetKEGG.PathLib(mSet, 'hsa', 'current')
mSet <- SetMetabolomeFilter(mSet, TRUE)             # TRUE alone does not restrict the background --
mSet <- Setup.KEGGReferenceMetabolome(mSet, 'reference_metabolome.txt')  # this call does (see pathway-mapping)
mSet <- CalculateOraScore(mSet, 'rbc', 'hyperg')
# Sends the mapped compound list to xialab.ca and can reject a filtered request outright on this
# version; if so, use pathway-mapping's Local-Only ORA (KEGGREST + local phyper, no remote call).

# Path B: no IDs -> mummichog on the FULL peak table (m/z + p-value + t-score)
# mSet <- InitDataObjects('mass_all', 'mummichog', FALSE)
# mSet <- UpdateInstrumentParameters(mSet, 5.0, 'negative')   # ppm + ionization mode are mandatory
# mSet <- Read.PeakListData(mSet, 'peaks_full.txt')           # ENTIRE table, not significant-only
# mSet <- PerformPSEA(mSet, 'hsa_mfn', 'current', permNum = 1000)
```

## Alternative Front End -- MS-DIAL

When peak detection happens in the MS-DIAL GUI/console (MS2Dec deconvolution, GC-EI, DIA/SWATH), import the alignment-result table and enter the pipeline at Stage 2. The framing is unchanged: the imported table is still a parameterized hypothesis. See metabolomics/msdial-preprocessing for the export-parsing details, then continue with normalization-qc onward.

## QC Checkpoints

Each gate hands off to its component skill when it fails; "refresh" means re-run the upstream stage with revised parameters, not patch the symptom downstream.

| Stage | Keep (pass) | Refresh (fail) -> where |
|---|---|---|
| Feature extraction | EIC + alignment of top hits inspect cleanly; feature count plausible after redundancy collapse | Tune centWave/bw against EIC FWHM -> xcms-preprocessing |
| Drift correction | QC RSD dropped AND biological-sample RSD unchanged; dilution-QC linearity holds | Back off spline span / exclude weak-in-QC features -> normalization-qc |
| QC quality | QC RSD <= 20-30%, D-ratio <= 0.5, blank ratio >= 3-5x | Drop failing features; check instrument/injection -> normalization-qc |
| Missingness | wholly-NA samples (QCRSC below-minQC wipeout) dropped, not imputed; remaining sparse holes imputed by mechanism; zero NAs remain in the matrix Stage 4 receives | Detection-rate filter before imputing; drop wiped samples then log2/2^x round-trip QRILC -> normalization-qc |
| PCA / QC clustering | pooled QCs cluster tightly at center; no batch-driven separation | Revisit batch correction / design -> normalization-qc, experimental-design/batch-design |
| Multivariate | Q2 high AND pQ2 small (permI >= 1000); `getSummaryDF()` has >0 rows (not a silent empty OPLS-DA model); PCA shows the same structure | Do not report a noise-separated score plot; use `fit_discriminant_guarded` -> statistical-analysis |
| Annotation | each reported name carries an MSI level; ion families collapsed; every feature carries a `mode` column | Downgrade Level 3-5 names; do not promote a DB hit; never mix modes -> metabolite-annotation |
| Identified-ORA gate | `identified_compounds` built ONLY from Level 1-2 rows, checked in code, not just narrated | Filter by `msi_level` before `Setup.MapData` -> metabolite-annotation, pathway-mapping |
| Pathway background | ORA uses assay-coverage background; mummichog uses the FULL table | `SetMetabolomeFilter(TRUE)` + `Setup.KEGGReferenceMetabolome` / supply R_all; Local-Only ORA if the remote call is rejected -> pathway-mapping |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `could not find function "readMSData"` | Legacy xcms <3 API | Use `readMsExperiment()` + `*Param` verbs (xcms 4.x) |
| `unused argument (ppm = ...)` | Loose args to `findChromPeaks` | Wrap in `CentWaveParam(...)`, pass via `param =` |
| Features on uncorrected RT | Grouped before alignment and never grouped after | Group AFTER `adjustRtime` (obiwarp needs no pre-grouping; PeakGroups alignment needs group -> align -> regroup) |
| Significance explodes after imputation | Half-min impute then test | Mechanism-aware QRILC/GSimp on sparse holes only |
| Clean OPLS-DA plot but it is noise | `permI = 20` (ropls default) | `permI >= 1000`; read `pQ2` from `getSummaryDF` |
| FDR is actually Holm/Holm-Sidak | R `p.adjust` default `'holm'`; statsmodels `'hs'` | Pass BH explicitly |
| Every pathway is significant | ORA on all-of-KEGG / mummichog on significant-only | Assay-coverage background; supply the FULL feature table |
| `opls() ERROR: missing value where TRUE/FALSE needed` at Stage 4, on real multi-batch data | Stage 2's imputation was skipped (a comment, not a call), OR `QCRSC` wiped an entire batch (fewer than `minQC` QCs) to all-NA and it reached `opls()` unimputed | Run Stage 2's drop-wiped-samples + `impute.QRILC` steps above and assert `!anyNA(imputed)` before Stage 4; never hand `normalized` (pre-drop, pre-imputation) to `opls()` |
| `impute.QRILC` itself errors `0 (non-NA) cases` / `NaNs produced` in `qnorm` | A sample (column) is 100% NA -- QRILC fits a per-sample quantile model and has nothing to fit when the whole column is missing | This is the QCRSC-wipeout case above, not a QRILC bug: drop the wholly-NA samples first, then QRILC on what remains (sparse holes only) |
| pmp step silently "works" after `fm <- t(feat)` but on features that look transposed | `featureValues()` already returns features-in-rows; the extra transpose inverts it and only survives because `pmp`'s `check_peak_matrix` auto-corrects (no warning when it can) | Skip the transpose (`fm <- feat`); add the `stopifnot(nrow(fm) == nrow(defs), ...)` check above so a real mismatch stops the pipeline instead of hiding |
| `identified_compounds` (or the mummichog fallback) is empty | Every annotated feature is Level 3-5; the MSI gate filtered all of them out, correctly | Do not loosen the filter -- report the honest coverage, or use Path B (mummichog on the full m/z table) instead |

## References

- Smith CA, Want EJ, O'Maille G, Abagyan R, Siuzdak G. 2006. XCMS: processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching, and identification. *Anal Chem* 78:779-787.
- Broadhurst D, Goodacre R, Reinke SN, Kuligowski J, Wilson ID, Lewis MR, Dunn WB. 2018. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. *Metabolomics* 14:72.
- Westerhuis JA, Hoefsloot HCJ, Smit S, Vis DJ, Smilde AK, et al. 2008. Assessment of PLSDA cross validation. *Metabolomics* 4:81-89.
- Schymanski EL, Jeon J, Gulde R, Fenner K, Ruff M, Singer HP, Hollender J. 2014. Identifying small molecules via high resolution mass spectrometry: communicating confidence. *Environ Sci Technol* 48:2097-2098.
- Wieder C, Frainay C, Poupin N, Rodriguez-Mier P, Vinson F, Cooke J, Lai RPJ, Bundy JG, Jourdan F, Ebbels T. 2021. Pathway analysis in metabolomics: recommendations for the use of over-representation analysis. *PLOS Comput Biol* 17(9):e1009105.

## Related Skills

- metabolomics/xcms-preprocessing - Stage 1 feature extraction parameters and the feature-table-as-artifact framing
- metabolomics/normalization-qc - Stage 2 drift correction, RSD/D-ratio filtering, PQN, mechanism-aware imputation
- metabolomics/metabolite-annotation - Stage 3 MSI/Schymanski confidence levels
- metabolomics/statistical-analysis - Stage 4 permutation-validated multivariate and dependence-aware FDR
- metabolomics/pathway-mapping - Stage 5 ORA vs mummichog and background construction
- metabolomics/msdial-preprocessing - Alternative front end entering at Stage 2
- metabolomics/lipidomics - Lipid-specific peak widths and annotation
- metabolomics/targeted-analysis - Absolute quantification branch
- metabolomics/isotope-tracing - Separate stable-isotope flux branch, not a stage of this untargeted pipeline
- multi-omics-integration/mofa-integration - Integrating the feature table with other omics layers
