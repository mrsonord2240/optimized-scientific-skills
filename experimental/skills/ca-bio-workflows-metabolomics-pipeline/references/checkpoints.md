## QC Checkpoints

Each gate hands off to its component skill when it fails; "refresh" means re-run the upstream stage with revised parameters, not patch the symptom downstream.

| Stage | Keep (pass) | Refresh (fail) -> where |
|---|---|---|
| Feature extraction | EIC + alignment of top hits inspect cleanly; feature count plausible after redundancy collapse | Tune centWave/bw against EIC FWHM -> xcms-preprocessing |
| Drift correction | QC RSD dropped AND biological-sample RSD unchanged; dilution-QC linearity holds | Back off spline span / exclude weak-in-QC features -> normalization-qc |
| QC quality | QC RSD <= 20-30%, D-ratio <= 0.5, blank ratio >= 3-5x | Drop failing features; check instrument/injection -> normalization-qc |
| Missingness | wholly-NA samples (QCRSC below-minQC wipeout) dropped, not imputed; remaining sparse holes imputed by mechanism; zero NAs remain in the matrix Stage 4 receives | Detection-rate filter before imputing; drop wiped samples then log2/2^x round-trip QRILC -> normalization-qc |
| PCA / QC clustering | pooled QCs cluster tightly at center; no batch-driven separation | Revisit batch correction / design -> normalization-qc, bio-experimental-design-batch-design |
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
