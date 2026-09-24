# Stage 1 -- Feature Extraction (modern xcms 4.x)

**Goal:** Turn centroided mzML into a features-by-samples table, carrying the parameters as part of the result.

**Approach:** Use the `MsExperiment`/`XcmsExperiment` containers with `*Param` objects; align to pooled QC, group AFTER alignment (obiwarp aligns the raw profile directly, so no pre-grouping is needed; the PeakGroups method instead needs group -> align -> regroup because it uses grouped anchor peaks), and treat filled values as imputations. Full parameter rationale (ppm, peakwidth, bw, prefilter) lives in metabolomics/xcms-preprocessing.

```r
# after defining mzml_files, pd, ionization_mode:
source('scripts/stage1_xcms_extract.R')   # -> xdata, feat (features x samples), defs
```
