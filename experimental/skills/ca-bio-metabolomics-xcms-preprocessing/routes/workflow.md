# File-based preprocessing

```bash
Rscript scripts/process_files.R sample_metadata.csv output-directory
```

CSV requires ordered `file_path,sample_name,sample_group` columns for centroided LC-MS files. The wrapper adapts the source example's broad-HPLC settings: centWave ppm 25, peakwidth 20–80 seconds, snthresh 10, prefilter 3/1000, noise 1000; obiwarp binSize 0.6; grouping bw 30, minFraction 0.5, binSize 0.025. These are example settings, not instrument-independent defaults; inspect EICs and adapt before use. It writes unfilled and filled intensity matrices plus a saved XCMS object and feature definitions. Existing names may be overwritten. It does not generate a per-feature filled-fraction report; derive it from the saved object's peak associations and flags before inference. Residual missing values remain possible.

Done when sample order, detection, alignment, grouping, fill provenance and output parameters are inspected and recorded.
