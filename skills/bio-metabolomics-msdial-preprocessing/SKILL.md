---
name: bio-metabolomics-msdial-preprocessing
description: Runs the MS-DIAL preprocessing workflow (peak picking, MS2Dec spectral deconvolution, alignment, gap-filling) and imports the alignment-result table into R or Python with honest filtering. Use when preprocessing LC-MS DDA/DIA (SWATH) or GC-MS raw data with MS-DIAL, deciding MS-DIAL vs XCMS, configuring the MSDIALCUI console run, or parsing an MS-DIAL export into a clean feature matrix. For programmatic R peak detection and the feature-table-as-artifact framing see metabolomics/xcms-preprocessing; for lipid annotation mode see metabolomics/lipidomics; for MSI-level confidence honesty see metabolomics/metabolite-annotation; for drift correction and QC see metabolomics/normalization-qc.
tool_type: mixed
primary_tool: msdial
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: MS-DIAL console 5.5.260820 (`MSDIALCUI.exe` -- confirmed installed and run on real data 2026-09-16), pandas 2.2+, R 4.3+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: run `MSDIALCUI.exe` (or `MSDIALCUI` on Linux/Mac) with no arguments to print the current subcommand/flag list, then `MSDIALCUI.exe <subcommand> --help` for that subcommand's real options
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show pandas` then `help(module.function)` to check signatures

The MS-DIAL GUI (`MSDIAL.exe`) runs only on Windows and hangs waiting for input if launched with `--help` or any unrecognized argument; the **separately-distributed** console app (`MSDIALCUI.exe`, packaged as its own `MSDIAL.console.*` release asset) is the cross-platform headless entry, and its real subcommand list -- confirmed by running it with no arguments -- is `gcms`, `lcms`, `lcimms`, `dims`, `imms`, `msn`, `eic`, `rtcorrection`, `imagegen`. GC-MS is real and working here: `gcms` is backed by its own `MsdialGcMsApi.dll` in the same release, so GC-EI does not require sourcing a separate MS-DIAL 4 build (see "Why GC-EI Is Different" below for what changed). If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt rather than retrying.

# MS-DIAL Preprocessing

**"Process my LC-MS run with MS-DIAL and give me a feature table"** -> Pick peaks per file, deconvolve chimeric MS/MS into clean component spectra (MS2Dec), align across samples, gap-fill, then import the alignment result and filter it honestly.
- CLI: `MSDIALCUI.exe lcms -i <in> -o <out> -m <param.txt>` for LC-MS (DDA and DIA share this one subcommand -- see "Run MS-DIAL Headless" for how acquisition type is actually set); `MSDIALCUI.exe gcms -i <in> -o <out> -m <param.txt>` for GC-MS
- R: `read.csv(..., skip = 4, check.names = FALSE)` to parse the alignment export
- Python: `pandas.read_csv(..., skiprows=4)` for the same export

## The Single Most Important Insight -- Preprocessing Software Is Not Neutral

The same raw files through MS-DIAL versus XCMS yield different feature tables and different marker lists. Li 2018 benchmarked five tools on a 1,100-compound standard and found that while feature *detection* was broadly similar, *quantification* and the set of selected discriminating markers differed by tool. A metabolomics "hit" is conditional on (raw data + software + version + every parameter + fill/filter order), not on the raw files alone. MS-DIAL's specific differentiator is **MS2Dec deconvolution**: it reconstructs clean, library-matchable MS/MS spectra from chimeric DDA/DIA fragment data, which is what makes wide-window DIA (SWATH) tractable at all. Report the full processing specification as part of the result, and treat a finding that survives only one pipeline as a candidate, not a result.

## MS-DIAL vs XCMS

| Axis | MS-DIAL | XCMS |
|---|---|---|
| Interface | Windows GUI + cross-platform console | R package (scriptable everywhere) |
| Core differentiator | MS2Dec MS/MS deconvolution (DDA + DIA) | centWave peak picking, full programmatic control |
| Annotation | Built-in (library + MS-FINDER + LipidBlast) | Separate (CAMERA, downstream tools) |
| Lipidomics | Strong (predicted-CCS / EAD structural elucidation in v5) | Manual |
| Reproducibility unit | Param file + GUI choices | Versioned R script |
| Best when | DIA data, lipidomics, GUI workflow, built-in IDs | Scripted pipelines, custom parameters, cohort scale |

Use MS-DIAL when DIA deconvolution or built-in lipid annotation is the point; use metabolomics/xcms-preprocessing for fully scripted, version-pinned cohort processing. The strongest untargeted claims replicate across both.

## Decision Tree by Scenario

| Situation | Do | Why |
|---|---|---|
| LC-MS, top-N MS/MS (DDA) | `MSDIALCUI.exe lcms` (folder input, or `acquisition_type=DDA` per file in a CSV `-i`) | Cleaner per-precursor MS2, but intensity-biased, stochastic coverage |
| LC-MS, wide-window MS/MS (DIA / SWATH) | `MSDIALCUI.exe lcms` (ABF input only; set `acquisition_type=DIA` per file in a CSV `-i`) | Complete MS2 coverage; chimeric spectra REQUIRE MS2Dec to be usable |
| GC-EI run | `MSDIALCUI.exe gcms`, or AMDIS/eRah | EI fragments every co-eluting compound; deconvolution IS detection (see below) |
| Headless / Linux cluster | `MSDIALCUI` (cross-platform console binary) with a `-m` param file | GUI is Windows-only; console is the reproducible batch path |
| Lipid-focused study | MS-DIAL + LipidBlast | -> metabolomics/lipidomics for lipid annotation mode |
| Fixed panel of known targets; MRM/SRM/PRM quantification | Do not use this untargeted-alignment workflow | -> metabolomics/targeted-analysis for calibration, internal standards, and validated concentrations; the console's `-t`/`--target` mode is not a substitute for a targeted assay |
| Already have an alignment CSV | skip processing, parse + filter | See import + honest-filter sections below |

## Why GC-EI Is Different

In GC-EI, 70 eV ionization fragments every compound reproducibly, so the trace at any retention time is a superposition of fragments from several co-eluting molecules. Naive peak picking conflates them; **deconvolution into component spectra IS the feature-detection step**, then each component is matched against EI+RI libraries (NIST, FiehnLib). Cross-run/cross-lab alignment uses **retention index** (Kovats n-alkanes, or Fiehn FAME markers giving diagnostic m/z 74/87) rather than raw RT, because RT drifts with column aging. GC-MS is supported directly by the current MS-DIAL 5.x console's `gcms` subcommand (confirmed real and working on 5.5.260820, backed by its own `MsdialGcMsApi.dll`) -- a separate MS-DIAL 4 build is not required.

## Run MS-DIAL Headless (console)

**Goal:** Process a folder of converted spectra into an alignment table without the GUI.

**Approach:** Pick the real subcommand (`lcms` or `gcms` -- there is no separate `lcmsdda`/`lcmsdia` token in the current console), point `-i`/`-o`/`-m` at input dir, output dir, and a method (parameter) file; keep `-p` only if the project should reopen in the GUI.

Verified 2026-09-16 against the real, installed `MSDIALCUI.exe` 5.5.260820: `MSDIALCUI.exe lcms -i <dir> -o <dir> -m <method.txt>` ran end to end on a real DDA mzML file and produced a real alignment export (see "Import the Alignment Result" below for its actual filename and column layout).

```bash
# LC-MS (DDA or DIA -- one subcommand covers both). Accepts netCDF/mzML/ABF for DDA;
# DIA/SWATH accepts ABF ONLY (convert vendor raw -> ABF first; MS2Dec is the point).
MSDIALCUI.exe lcms -i ./LCMS/ -o ./LCMS_out/ -m ./Msdial-lcms-Param.txt

# Acquisition type (DDA vs DIA) is not a console token -- it is set per file. For a
# single-mode batch, a plain input folder is enough (every file processed the same
# way, as above). To mix modes or set it explicitly, pass -i a CSV instead of a folder,
# with one row per raw file and an `acquisition_type` column (DDA or DIA), per the
# official MS-DIAL 5 console tutorial (systemsomicslab.github.io/msdial5tutorial):
#   file_path,file_name,file_type,class_id,acquisition_type,batch_order,analytical_order,factor
#   LCMS/sample1.mzML,sample1,Sample,A,DDA,1,1,1
MSDIALCUI.exe lcms -i ./filelist.csv -o ./LCMS_out/ -m ./Msdial-lcms-Param.txt

# GC-EI: retention-index alignment, quant-mass quantification. Real and working in
# the current 5.x console (see "Why GC-EI Is Different" above) -- no MS-DIAL 4 needed.
MSDIALCUI.exe gcms -i ./GCMS/ -o ./GCMS_out/ -m ./Msdial-GCMS-Param.txt -p
```

The per-file `acquisition_type` mechanism above is documented by the official console tutorial and the CSV input itself has been accepted by the console, but it has not yet been exercised here on real DIA/ABF data. Treat true DIA MS2Dec behavior as a run-specific acceptance check: retain the input CSV, method file, and representative deconvolved spectra, then confirm them before interpreting DIA identifications.

The parameter (`-m`) file is plain text with **`Key: Value`** pairs (colon-separated, one per line, `#` for comments) -- not `Key=Value`. A malformed or misspelled line can be **silently ignored**: the console may exit 0 with no warning while using its built-in default. Before trusting a batch, run a small representative input twice with one safely extreme, known key (for example `Minimum peak height`) and confirm that the output feature count changes in the expected direction; archive both parameter files and counts. The `Minimum peak height` key is the direct analog of an intensity floor and is instrument-dependent: the GUI default is tuned for a TOF and is often far too high (or its baseline assumption wrong) for an Orbitrap. Set the alignment reference to a pooled QC, never to file #1 by default.

## Import the Alignment Result into R

**Goal:** Split the MS-DIAL alignment export into a feature-metadata frame and an intensity matrix.

**Approach:** The export carries four header rows above the real column header (sample class / file type / injection order / batch), so skip them. Verified 2026-09-16 against a real `MSDIALCUI.exe` 5.5.260820 run: the file is **not** named `AlignResult.txt` -- the console writes `AlignResult-<timestamp>.mdalign` (still a plain tab-separated text file, plus sibling `.mzTab` and `.qa.tsv` exports of the same run); glob for `AlignResult-*.mdalign` rather than assuming the literal name. Metadata columns precede the per-sample intensity columns, but real exports interleave extra annotation/QC columns (Formula, Ontology, INCHIKEY, SMILES, matched-flag and dot-product-score columns, spectra text, ...) between the last *known* metadata column and the real sample columns -- taking "everything after the last known metadata column" pulls those text columns into what should be a numeric matrix. Anchor on the header block instead: MS-DIAL writes each of the 4 header rows as one row-label cell (`Class`, `File type`, `Injection order`, `Batch ID`) immediately followed by one value cell per sample, so the column right after the `Class` cell in row 1 is where the real sample columns start.

```r
# MS-DIAL alignment export: real column header is on row 5, so skip the first 4 rows.
# Filename is timestamped -- Sys.glob() instead of a literal name.
export_file <- Sys.glob('AlignResult-*.mdalign')[1]
msdial <- read.csv(export_file, sep = '\t', skip = 4, check.names = FALSE)

# Metadata columns appear before the per-sample intensity columns. Common ones:
# 'Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name', 'Adduct type',
# 'Fill %', 'MS/MS assigned', 'Reference RT', 'Formula', 'Ontology', 'INCHIKEY',
# 'SMILES', 'Annotation tag (VS1.0)'.
meta_cols <- c('Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name',
               'Adduct type', 'Fill %', 'MS/MS assigned', 'Annotation tag (VS1.0)')
meta_cols <- intersect(meta_cols, colnames(msdial))

# Sample columns start right after the header block's "Class" label cell -- NOT
# simply "everything after the last known metadata column" (real exports interleave
# unlisted annotation/QC columns before the real per-sample data; see above).
header_row1 <- as.character(read.table(export_file, sep = '\t', header = FALSE, nrows = 1,
                                        check.names = FALSE, colClasses = 'character',
                                        comment.char = '')[1, ])
class_idx <- which(header_row1 == 'Class')
sample_cols <- colnames(msdial)[(class_idx + 1):ncol(msdial)]

feature_info <- msdial[, meta_cols]
intensity <- as.matrix(msdial[, sample_cols])
rownames(intensity) <- msdial[['Alignment ID']]
storage.mode(intensity) <- 'numeric'  # confirms every detected sample col is truly numeric
```

## Import the Alignment Result into Python

**Goal:** Same split, in pandas.

**Approach:** `skiprows=4` to land on the real header; locate the sample columns via the same `Class`-cell anchor as the R path above, not by position after the last known metadata column.

```python
import csv
import pandas as pd
import glob

export_file = glob.glob('AlignResult-*.mdalign')[0]  # timestamped filename, not a literal AlignResult.txt

with open(export_file, newline='') as f:
    header_row1 = next(csv.reader(f, delimiter='\t'))
class_idx = header_row1.index('Class')

msdial = pd.read_csv(export_file, sep='\t', skiprows=4)
meta_cols = ['Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name', 'Adduct type', 'Fill %', 'MS/MS assigned', 'Annotation tag (VS1.0)']
meta_cols = [c for c in meta_cols if c in msdial.columns]
sample_cols = msdial.columns[class_idx + 1:]  # everything after the "Class" cell, not after the last known metadata column

feature_info = msdial[meta_cols].copy()
intensity = msdial[sample_cols].astype(float).copy()
intensity.index = msdial['Alignment ID']
```

## Filter the Table Honestly

**Goal:** Keep features supported by real signal and known confidence, without overtrusting annotation tags.

**Approach:** Filter on Fill% (cross-sample presence), require MS/MS support for any feature called identified, and tie the annotation tag to a real MSI confidence level rather than treating a name as proof.

```r
# Fill% is the fraction of samples with a DETECTED (not gap-filled) peak. Low Fill% means
# the feature exists mostly as gap-filled noise-floor integrals, which fabricate intensity
# (an honest 'below detection' becomes a positive number). 70% is a common floor -- but
# MS-DIAL 5.x's real "Fill %" column is a 0-1 FRACTION (e.g. 1.00 = 100%), not 0-100
# (confirmed on a real MSDIALCUI.exe 5.5.260820 run, 2026-09-16), so compare to 0.70:
keep_fill <- feature_info[['Fill %']] >= 0.70

# An annotated name without MS/MS is at best an MSI Level 2/3 putative ID (accurate mass
# only). Require MS/MS support before trusting any identity downstream. MS-DIAL 5.x's real
# value is title-case text ("True"/"False"), not "TRUE"/"FALSE" -- an exact-case comparison
# silently keeps zero features on real output, so compare case-insensitively:
has_msms <- tolower(trimws(feature_info[['MS/MS assigned']])) == 'true'

# Annotation tag confidence (do NOT treat a name as an identification). The exact tag
# vocabulary is MS-DIAL-version-dependent, so inspect unique(feature_info[['Annotation tag (VS1.0)']])
# and map the strings the build actually emits rather than hard-coding them:
#   Metabolite / Lipid  with MS/MS  -> MSI Level 2 (spectral library match)
#   Suggested*          mass-only   -> MSI Level 3 (putative, no MS/MS)
#   Unknown                         -> unannotated feature
feature_info$msi_level <- ifelse(feature_info[['Annotation tag (VS1.0)']] %in% c('Metabolite', 'Lipid') & has_msms, 2,
                          ifelse(grepl('^Suggested', feature_info[['Annotation tag (VS1.0)']]), 3, NA))

filtered <- intensity[keep_fill, ]
```

**Python:** the same `MS/MS assigned` check needs a different guard, not just a different case. `pandas.read_csv` auto-casts a text column of literal `True`/`False` (or `TRUE`/`FALSE`) to native `bool` with no warning, so the R idiom's natural Python port (`== 'True'`) silently evaluates to all-`False` -- confirmed 2026-09-16 (0/94 real MS/MS-supported features flagged when ported naively). Force the comparison to be dtype-safe:

```python
keep_fill = feature_info['Fill %'] >= 0.70  # same 0-1 fraction as the R path
has_msms = feature_info['MS/MS assigned'].astype(str).str.strip().str.lower() == 'true'
```

Confidence-level honesty and orthogonal-evidence identification belong to metabolomics/metabolite-annotation; this skill only routes the tag to the right level. Fill% / blank / drift filtering interacts with normalization-qc - process blanks and pooled QCs through the SAME run, then filter the aligned table.

## Per-Method Failure Modes

### DIA processed as DDA (wrong acquisition type)
- **Trigger:** Running SWATH/DIA data through `lcms` without ever setting DIA acquisition type.
- **Mechanism:** Without DIA acquisition type set (per file, via a CSV `-i` with an `acquisition_type` column -- see "Run MS-DIAL Headless"), MS-DIAL does not deconvolve wide-isolation chimeric MS/MS, so fragments from co-isolated precursors stay mixed.
- **Symptom:** Library matches to the wrong compound; "clean-looking" spectra that fail orthogonal confirmation.
- **Fix:** Pass `-i` a CSV with `acquisition_type=DIA` for those files (ABF input only); MS2Dec deconvolution is the entire reason to run DIA in MS-DIAL.

### Over-trusting the annotation tag
- **Trigger:** Filtering on `Annotation tag != Unknown` and calling the survivors "identified."
- **Mechanism:** A `Suggested*` tag is an accurate-mass guess with no MS/MS; a named hit without MS/MS is MSI Level 3.
- **Symptom:** A marker list full of confident-sounding names that do not validate against standards.
- **Fix:** Require MS/MS support for any identity claim (case-insensitive check -- MS-DIAL 5.x's real value is `True`/`False`, not `TRUE`/`FALSE`); map tags to MSI levels (see filtering section) and defer to metabolomics/metabolite-annotation.

### Gap-fill masquerading as measurement
- **Trigger:** Treating low-Fill% features as quantitative.
- **Mechanism:** Gap-filling integrates whatever signal sits in the m/z-RT box even when no peak exists, turning a true below-detection (MNAR, left-censored) value into a positive number.
- **Symptom:** "Significant" features that are mostly gap-filled in one group; shrunken fold-changes for on/off markers.
- **Fix:** Report per-feature filled fraction; gate on Fill% >= 0.70 (MS-DIAL 5.x reports it as a 0-1 fraction, not 0-100 -- see filtering section); for inferential stats prefer MNAR-aware imputation over naive fill (see metabolomics/normalization-qc).

### GC-EI run through an LC pipeline
- **Trigger:** Treating GC-EI data like LC peak-pick-then-group, or assuming MS-DIAL 5 has no GC mode.
- **Mechanism:** EI needs component deconvolution, not adduct-style peak picking, and RI (not RT) alignment; the current MS-DIAL 5.x console's `gcms` subcommand does handle this (confirmed real and working on 5.5.260820) if a plain LC-style `lcms` run is used instead.
- **Symptom:** Conflated co-eluting compounds and cross-lab RT misalignment.
- **Fix:** Use `MSDIALCUI.exe gcms` (or AMDIS/eRah); align on Kovats/FAME retention index.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| Fill% >= 70% (0.70 on MS-DIAL 5.x's real 0-1 Fill % scale) | Common untargeted practice | Below this, the feature is mostly gap-filled noise-floor integrals, not measurements |
| QC CV (RSD) < 20-30% | Broadhurst 2018 | Technical reproducibility floor; drop features noisier than this in pooled QCs |
| D-ratio (sd_QC/sd_sample) < 0.5 | Broadhurst 2018 | Keeps features whose technical variance is well below biological variance |
| Blank filter: sample mean > 3-5x blank mean | Broadhurst 2018 | Removes background/contaminant features present in process blanks |
| ~10x more features than compounds | Mahieu 2017 | One metabolite makes adducts/isotopes/fragments; counting features over-counts hypotheses |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `MsdialConsoleApp lcmsdda`/`lcmsdia` not found | That naming is from the older pre-5.x console; current builds ship `MSDIALCUI.exe` with one `lcms` subcommand | Run `MSDIALCUI.exe` with no args to print the real subcommand list; use `lcms` for both DDA and DIA |
| All columns land in one field on import | Header offset wrong; tab-separated export read as CSV | `skip=4` (R) / `skiprows=4` (Python), set `sep='\t'` |
| `Fill %` filter keeps 0 features, or keeps everything | Comparing a 0-1 fraction (MS-DIAL 5.x's real scale) against a 70-100 threshold | Compare against 0.70, not 70 |
| `MS/MS assigned` filter always evaluates False | Real value is title-case `True`/`False` text (Python: pandas may also silently cast it to native `bool`), not `TRUE`/`FALSE` | Compare case-insensitively: `tolower(trimws(x)) == 'true'` (R) / `x.astype(str).str.strip().str.lower() == 'true'` (Python) |
| DIA input rejects mzML | DIA mode accepts ABF only | Convert vendor raw to ABF (Reifycs ABF converter) before setting `acquisition_type=DIA` |
| `Annotation tag` column not found | Header changes across versions (e.g. `Annotation tag (VS1.0)`) | Match by prefix / inspect `colnames()`; do not hard-code the suffix |
| Console command not found on Linux | Expecting the GUI executable | The GUI (`MSDIAL.exe`) is Windows-only and hangs on unrecognized args; run the separately-packaged `MSDIALCUI` console binary (cross-platform) |
| Few features detected | `Minimum peak height` default too high for the instrument | Lower it toward the real baseline; defaults are TOF-tuned |
| Clean console exit but a method change appears to have no effect | A `Key=Value` line or misspelled method key was silently ignored, so defaults were used | Use only `Key: Value`; prove the method is read by changing one known key on a representative input and checking that the resulting feature count changes, then retain both methods and counts |

## References

- Tsugawa H, Cajka T, Kind T, Ma Y, Higgins B, Ikeda K, Kanazawa M, VanderGheynst J, Fiehn O, Arita M. MS-DIAL: data-independent MS/MS deconvolution for comprehensive metabolome analysis. *Nat Methods.* 2015; 12(6):523-526.
- Tsugawa H, Ikeda K, Takahashi M, et al. A lipidome atlas in MS-DIAL 4. *Nat Biotechnol.* 2020; 38(10):1159-1163.
- Takeda H, Takahashi M, Ikeda K, et al. MS-DIAL 5 multimodal mass spectrometry data mining unveils lipidome complexities. *Nat Commun.* 2024; 15:9903.
- Li Z, Lu Y, Guo Y, Cao H, Wang Q, Shui W. Comprehensive evaluation of untargeted metabolomics data processing software in feature detection, quantification and discriminating marker selection. *Anal Chim Acta.* 2018; 1029:50-57.
- Mahieu NG, Patti GJ. Systems-level annotation of a metabolomics data set reduces 25,000 features to fewer than 1,000 unique metabolites. *Anal Chem.* 2017; 89(19):10397-10406.
- Broadhurst D, Goodacre R, Reinke SN, Kuligowski J, Wilson ID, Lewis MR, Dunn WB. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. *Metabolomics.* 2018; 14(6):72.
- Stein SE. An integrated method for spectrum extraction and compound identification from gas chromatography/mass spectrometry data (AMDIS). *J Am Soc Mass Spectrom.* 1999; 10(8):770-781.

## When NOT to Use This Skill

- Need fully scripted, version-pinned R processing with custom peak-picking parameters -> metabolomics/xcms-preprocessing.
- The data is already lipid-annotated or the question is lipid-class/species-level -> metabolomics/lipidomics.
- The question is about confidence-level honesty or orthogonal-evidence identification for a feature you already have -> metabolomics/metabolite-annotation.
- The question is drift correction, QC-CV/D-ratio filtering, blank filtering, or MNAR-aware imputation on an already-aligned table -> metabolomics/normalization-qc (this skill only routes the Fill%/MS-MS/annotation-tag columns to the right level; it does not implement those filters).
- Need concentrations for a predefined MRM/SRM/PRM panel of known analytes -> metabolomics/targeted-analysis. MS-DIAL console target mode (`-t`/`--target`) does not replace calibration, co-eluting internal standards, and method validation.

## Related Skills

- metabolomics/xcms-preprocessing - Programmatic R preprocessing and the feature-table-as-artifact framing
- metabolomics/lipidomics - Lipid annotation mode and LipidBlast workflows
- metabolomics/metabolite-annotation - MSI confidence levels and orthogonal-evidence identification
- metabolomics/normalization-qc - Drift correction, QC/CV/D-ratio filtering, MNAR-aware imputation
- metabolomics/targeted-analysis - Validated targeted MRM/SRM/PRM quantification for a predefined analyte panel
