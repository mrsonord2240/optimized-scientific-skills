---
name: bio-proteomics-data-import
description: Loads mass-spectrometry data into Python/R and strips the search engine's bookkeeping before any number is trusted -- removes decoys (REV__/Reverse), contaminants (CON__/Potential contaminant), Only-identified-by-site groups, and resolves semicolon razor/leading protein-ID ambiguity in MaxQuant proteinGroups.txt, DIA-NN report.parquet, and mzML/mzXML. Distinguishes Intensity (raw) vs LFQ intensity (MaxLFQ) vs iBAQ, treats a MaxQuant zero as missing (NaN, not log2(-inf)), and diagnoses the missingness contract (intensity-dependent in both DDA and DIA; DIA has fewer missing values). Use when starting an analysis from raw spectra or a search engine output. Downstream normalization and stats are differential-abundance; reporter-ion/MaxLFQ quant is quantification; protein grouping is protein-inference.
tool_type: mixed
primary_tool: pyOpenMS
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: pyOpenMS 3.5.0, pandas 3.0.5, numpy 2.5.3 (Python blocks, 2026-09-15); MSnbase 2.28+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Mass Spectrometry Data Import -- Inheriting the Acquisition Contract and Stripping the Bookkeeping

**"Load my mass spec data into Python"** -> Parse spectra or a search-engine table AND immediately enforce two contracts -- which quant column carries real biology, and which rows are search-engine bookkeeping that must be deleted -- because the same proteinGroups.txt yields different conclusions depending on the column read and the rows kept.
- Python: `pyopenms.MzMLFile().load(path, exp)` for `.mzML`, `MzXMLFile().load(path, exp)` for `.mzXML`; `pandas.read_csv(sep='\t')` for MaxQuant; `pandas.read_parquet` for DIA-NN
- R: `Spectra::Spectra()` / `QFeatures::readQFeatures()` for raw and quantified data (MSnbase still works but is in maintenance mode)

Scope: this skill owns reading spectra/search outputs into memory, deleting decoy/contaminant/site-only rows, selecting an already-computed quant column, and characterizing missingness. Format conversion (RAW -> mzML) -> peptide-identification. MaxLFQ/TMT reporter **computation or correction** -> quantification. Protein-group parsimony -> protein-inference. Normalization and imputation -> differential-abundance and expression-matrix/normalization. OUT OF SCOPE: statistical testing, batch correction, and the actual imputation step (this skill only diagnoses the missingness so the right imputer is chosen later).

## The Single Most Important Modern Insight -- Import Is Where Two Contracts Are Read and Enforced

1. **A "data import" is never just file parsing -- it is the moment the acquisition mode's quantitative contract and its missingness structure are inherited.** DDA selects the top-N most intense precursors per cycle, and which precursors get picked is partly stochastic and abundance-biased, so the same low-abundance peptide is sampled in run A and missed in run B; this manufactures structured, left-censored MNAR missingness. DIA fragments every precursor in every window every cycle, so it has fewer missing values, but those it has are still mostly intensity-dependent (concentrated in low-abundance proteins), not MCAR (Hediyeh-zadeh 2023). The catastrophic error this prevents: imputing a matrix with a mean/KNN method that assumes MCAR, which biases low-abundance proteins upward and manufactures false hits. The mode is born at acquisition and inherited at import, but the imputer is chosen from the missingness diagnosis made here, not from the acquisition mode alone.

2. **The search engine's bookkeeping must be stripped before any number is trusted.** A proteinGroups.txt carries decoy rows (`Reverse == '+'`, `REV__` prefix in the ID) from the target-decoy FDR machinery, contaminant rows (`Potential contaminant == '+'`, `CON__` prefix), and Only-identified-by-site rows (the protein has no unmodified-peptide evidence, only a modified site). Keeping any of these leaks non-biological signal into the intensity matrix and inflates IDs. The catastrophic error: reporting differential abundance on a matrix where decoy or keratin rows survived.

3. **The same proteinGroups.txt yields different biology from different columns, and a zero is not a measurement.** `Intensity` is raw summed precursor signal (not normalized, not comparable across samples for ratios). `LFQ intensity` is MaxLFQ-normalized and is the column for between-sample comparison. `iBAQ` is intensity divided by the number of observable tryptic peptides -- a within-sample molar proxy, not a between-sample quant. MaxQuant writes 0 for "not quantified", so log2(0) = -inf; replace 0 -> NaN before any transform. The catastrophic error: log2-transforming raw `Intensity` (or iBAQ) and reading the ratios as biology.

## Tool Taxonomy

| Tool / method | Citation | Mechanism / role | When |
|---|---|---|---|
| pyOpenMS `MzMLFile` / `MzXMLFile` | Rost 2014 | Loads the matching mzML/mzXML format into an MSExperiment; iterate spectra by MS level | Programmatic access to raw peaks, precursor m/z, isolation windows |
| pandas `read_csv`/`read_parquet` | -- | Tabular ingest of MaxQuant TSV and DIA-NN parquet | All search-engine output tables |
| DIA-NN report | Demichev 2020 | Long-format precursor table; `report.parquet` is the default (1.9+) and the only default (2.0) | DIA quant; pivot on `PG.MaxLFQ` after q-filtering |
| MaxQuant `txt/` outputs | Cox 2014 (MaxLFQ) | `proteinGroups.txt` (group level), `evidence.txt` (per-PSM) | DDA label-free; corrected TMT reporter columns can be imported but reporter computation routes out |
| Spectra + QFeatures (R) | -- | Current Bioconductor raw + quantified-feature containers; `readQFeatures`, `aggregateFeatures` | R pipelines; preferred over MSnbase going forward |
| MSnbase `readMSData` (R) | -- | On-disk raw reading; maintenance mode (route OUT to Spectra/QFeatures) | Legacy R code only |
| ThermoRawFileParser / msconvert | Hulstaert 2020 / Chambers 2012 | RAW -> mzML conversion (route OUT) | File conversion is peptide-identification |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|---|---|---|
| MaxQuant DDA label-free, between-sample comparison | Read `LFQ intensity` columns from proteinGroups.txt | MaxLFQ-normalized; the only MaxQuant column valid for cross-sample ratios |
| MaxQuant, absolute/molar abundance within one sample | Read `iBAQ` columns | iBAQ is a within-sample molar proxy; do not use across samples |
| Need raw uncorrected signal for a custom normalization | Read `Intensity` columns, normalize yourself | `Intensity` is raw summed precursor area, not comparable as-is |
| MaxQuant TMT proteinGroups.txt, already reporter-corrected | Read `Reporter intensity corrected <channel>` columns; map channels to samples from the experimental design | The lone `Intensity` column is a group total, not a per-channel matrix. Reporter correction/computation and channel-normalization choices -> quantification |
| DIA-NN output (1.9 or 2.0) | `pd.read_parquet('report.parquet')`, filter run-level and `Global.PG.Q.Value` q-values, pivot `PG.MaxLFQ`, 0 -> NaN, log2 | 2.0 dropped the TSV default; without the global protein q-value filter, groups that pass only within single runs leak into the cross-run matrix |
| Raw spectra, need peaks/precursor/isolation window | `MzMLFile` for `.mzML`, `MzXMLFile` for `.mzXML` | Programmatic peak and isolation-window access for QC and co-isolation reasoning |
| R-based MaxQuant proteinGroups import | `readQFeatures`, then `make.names()` before `filterFeatures`, `zeroIsNA`, `logTransform` | QFeatures makes row-data names syntactic; use the runnable block below. `aggregateFeatures` is for peptide-to-protein aggregation, not proteinGroups |
| Data came from DDA, planning imputation | Diagnose missingness as MNAR -> route to left-censored imputation | DDA top-N sampling makes missingness abundance-dependent |
| Data came from DIA, planning imputation | Run the same diagnostic; a negative abundance-missingness correlation means left-censored handling, as for DDA | DIA has fewer missing values, but they are still mostly intensity-dependent (Hediyeh-zadeh 2023) |

Default when uncertain: read `LFQ intensity` (label-free MaxQuant) or `PG.MaxLFQ` after q-filtering (DIA-NN), strip Reverse/contaminant/site-only rows, set 0 -> NaN, then diagnose missingness before choosing an imputer. For TMT, do not substitute the single `Intensity` total: use already-corrected reporter channels and route reporter processing to quantification.

## Loading mzML/mzXML with pyOpenMS

**Goal:** Parse raw spectra into memory for QC, peak access, and isolation-window reasoning.

**Approach:** Load into an MSExperiment (filled in place), iterate by MS level; `get_peaks()` returns a tuple of (mz, intensity) numpy arrays, and `getPrecursors()` returns a list that is empty for all-ion (AIF/MSE/bbCID) MS2 scans. An isolation width of 0 means the offsets were not written, not a 0-Th window.

```python
from pathlib import Path
from pyopenms import MSExperiment, MzMLFile, MzXMLFile

exp = MSExperiment()
path = 'sample.mzML'
suffix = Path(path).suffix.lower()
if suffix == '.mzml':
    loader = MzMLFile()
elif suffix == '.mzxml':
    loader = MzXMLFile()
else:
    raise ValueError('Expected an .mzML or .mzXML file')
loader.load(path, exp)  # fills exp in place; returns None

for spectrum in exp:
    if spectrum.getMSLevel() == 1:
        mz, intensity = spectrum.get_peaks()  # tuple of two numpy arrays
    elif spectrum.getMSLevel() == 2:
        precs = spectrum.getPrecursors()  # a list; empty for all-ion (AIF/MSE) MS2 scans
        if not precs:
            continue  # record as no-precursor MS2 instead of indexing [0]
        precursor = precs[0]
        precursor_mz = precursor.getMZ()
        window = precursor.getIsolationWindowLowerOffset() + precursor.getIsolationWindowUpperOffset()
        if window == 0:
            print(f'{spectrum.getNativeID()}: isolation offsets not written (width unknown, not 0 Th)')
```

## Loading and Cleaning MaxQuant proteinGroups.txt

**Goal:** Get a trustworthy log2 intensity matrix with bookkeeping rows removed and missing values represented as NaN.

**Approach:** Require the MaxQuant `Protein IDs` field, construct a row-indexed flag mask (so a table with all flag columns absent gets a useful error instead of `KeyError: True`), strip Reverse/contaminant/site-only rows, resolve semicolon IDs, pick `LFQ intensity` columns, set 0 -> NaN, and log2-transform. The counts are retained for the methods record.

```python
import pandas as pd
import numpy as np

pg = pd.read_csv('proteinGroups.txt', sep='\t', low_memory=False)  # mixed-type cols

if 'Protein IDs' not in pg:
    raise ValueError("Expected MaxQuant proteinGroups.txt with a 'Protein IDs' column; this is not that table")

# Keep a Series even if a bookkeeping field is absent. All three absent means this
# is not a MaxQuant proteinGroups table, so fail before selecting the quant family.
flag_cols = ('Reverse', 'Potential contaminant', 'Only identified by site')
if not any(c in pg for c in flag_cols):
    raise ValueError("Expected MaxQuant proteinGroups.txt bookkeeping columns; received none of Reverse, Potential contaminant, Only identified by site")
keep = pd.Series(True, index=pg.index)
for col in flag_cols:
    if col in pg:
        keep &= pg[col].fillna('').ne('+')
rows_read, rows_kept = len(pg), int(keep.sum())
pg = pg.loc[keep].copy()

# Protein IDs / Majority protein IDs / Gene names are SEMICOLON lists; take the first (leading/razor) entry
pg['leading_protein'] = pg['Protein IDs'].str.split(';').str[0]
pg['leading_gene'] = pg.get('Gene names', pd.Series('', index=pg.index)).fillna('').str.split(';').str[0]

lfq_cols = [c for c in pg.columns if c.startswith('LFQ intensity ')]  # MaxLFQ-normalized, between-sample comparable
if not lfq_cols:
    tmt_cols = [c for c in pg.columns if c.startswith('Reporter intensity corrected ')]
    if tmt_cols:
        raise ValueError("No LFQ intensity columns: this looks like MaxQuant TMT. Use 'Reporter intensity corrected <channel>' columns; reporter processing -> quantification")
    raise ValueError("No LFQ intensity columns: use a label-free MaxQuant proteinGroups.txt, or verify the supported quant family")
matrix = pg[['leading_protein', 'leading_gene'] + lfq_cols].copy()
matrix[lfq_cols] = matrix[lfq_cols].replace(0, np.nan)  # MaxQuant writes 0 for missing; log2(0) = -inf
matrix[lfq_cols] = np.log2(matrix[lfq_cols])
matrix = matrix[matrix[lfq_cols].notna().any(axis=1)]  # groups with no valid LFQ value carry no quant
print(f'MaxQuant rows: read={rows_read}, after bookkeeping={rows_kept}, quantified={len(matrix)}')
```

### Importing already-corrected MaxQuant TMT channels

Do not feed a TMT table to the LFQ block: its one `Intensity` field is a group total. This import uses only MaxQuant's already-corrected channel values; choose the reporter correction and channel normalization upstream in **quantification**, and supply the channel-to-sample mapping separately.

```python
import numpy as np
import pandas as pd

pg = pd.read_csv('proteinGroups.txt', sep='\t', low_memory=False)
flag_cols = ('Reverse', 'Potential contaminant', 'Only identified by site')
if 'Protein IDs' not in pg or not any(c in pg for c in flag_cols):
    raise ValueError('Expected MaxQuant proteinGroups.txt with Protein IDs and at least one bookkeeping column')
keep = pd.Series(True, index=pg.index)
for col in flag_cols:
    if col in pg:
        keep &= pg[col].fillna('').ne('+')
pg = pg.loc[keep].copy()
pg['leading_protein'] = pg['Protein IDs'].str.split(';').str[0]
pg['leading_gene'] = pg.get('Gene names', pd.Series('', index=pg.index)).fillna('').str.split(';').str[0]
reporter_cols = [c for c in pg.columns if c.startswith('Reporter intensity corrected ')]
if not reporter_cols:
    raise ValueError("No 'Reporter intensity corrected <channel>' columns; this is not a corrected MaxQuant TMT proteinGroups table")
tmt = pg[['leading_protein', 'leading_gene'] + reporter_cols].copy()
tmt[reporter_cols] = np.log2(tmt[reporter_cols].replace(0, np.nan))
tmt = tmt[tmt[reporter_cols].notna().any(axis=1)]
```

### MaxQuant proteinGroups.txt in R / QFeatures

`readQFeatures()` turns spaces in row-data columns into dots. Rename before filtering; `aggregateFeatures()` is not needed because the input is already protein-group level.

```r
library(QFeatures)
pg <- read.delim('proteinGroups.txt', quote = '', check.names = FALSE)
lfq <- grep('^LFQ intensity ', names(pg), value = TRUE)
stopifnot(length(lfq) > 0)
qf <- readQFeatures(pg, quantCols = lfq, name = 'proteinGroups')
rd <- rowData(qf[[1]]); colnames(rd) <- make.names(colnames(rd)); rowData(qf[[1]]) <- rd
qf <- filterFeatures(qf, ~ !(Reverse %in% '+') & !(Potential.contaminant %in% '+') & !(Only.identified.by.site %in% '+'))
qf <- zeroIsNA(qf, 1); qf <- logTransform(qf, i = 1, name = 'log2LFQ')
```

## Loading DIA-NN report.parquet

**Goal:** Reshape the long DIA-NN report into a confident protein-by-run matrix.

**Approach:** Read the parquet (default since 1.9, only default in 2.0), filter run-level precursor and protein-group q-values AND the experiment-wide `Global.PG.Q.Value` to 1% FDR BEFORE pivoting on `PG.MaxLFQ` (with library-based MBR on 1.9.x use `Lib.PG.Q.Value`), then set 0 -> NaN and log2.

```python
import numpy as np
import pandas as pd

report = pd.read_parquet('report.parquet')  # report.tsv dropped as default in DIA-NN 2.0
report = report[(report['Q.Value'] <= 0.01) & (report['PG.Q.Value'] <= 0.01)
                & (report['Global.PG.Q.Value'] <= 0.01)]  # run-level AND experiment-wide 1% FDR before quant

matrix = report.pivot_table(index='Protein.Group', columns='Run', values='PG.MaxLFQ', aggfunc='first')
matrix = np.log2(matrix.replace(0, np.nan))  # PG.MaxLFQ can be 0 too; log2(0) = -inf
```

## Diagnosing the Missingness Contract

**Goal:** Quantify the missing-value pattern so the legitimate imputation class can be chosen downstream.

**Approach:** Count NaN per protein and per sample on a log2 matrix with NaN for missing (zeros left in, or a linear scale, hide the signature). A negative correlation between missingness and mean abundance is the MNAR (left-censored) signature; expect it in DDA and, with fewer missing values, in DIA too.

```python
import numpy as np

def assess_missingness(matrix, sample_cols):
    miss_per_protein = matrix[sample_cols].isna().sum(axis=1)
    miss_per_sample = matrix[sample_cols].isna().sum(axis=0)
    total_pct = 100 * matrix[sample_cols].isna().sum().sum() / matrix[sample_cols].size
    mean_abund = matrix[sample_cols].mean(axis=1)  # negative corr with missingness => MNAR / left-censored
    mnar_corr = mean_abund.corr(miss_per_protein)
    return {'per_protein': miss_per_protein, 'per_sample': miss_per_sample, 'total_pct': total_pct, 'abundance_missing_corr': mnar_corr}
```

## Per-Method Failure Modes

### MaxQuant wrong quant column

**Trigger:** Reading `Intensity` (raw) or `iBAQ` when between-sample ratios are intended.
**Mechanism:** `Intensity` is un-normalized summed precursor signal; `iBAQ` is a within-sample molar proxy. Neither is comparable across samples the way `LFQ intensity` is.
**Symptom:** Ratios track total loaded protein / sample depth rather than biology; fold changes shift when one sample's loading changes.
**Fix:** Use `LFQ intensity` for cross-sample comparison; if computing custom normalization use `Intensity` and normalize explicitly (expression-matrix/normalization).

### TMT table sent to the LFQ route

**Trigger:** Passing a MaxQuant TMT `proteinGroups.txt` to code that expects `LFQ intensity <sample>` columns.
**Mechanism:** TMT stores per-channel values in `Reporter intensity corrected <channel>`; its lone `Intensity` is a group total and cannot replace a sample matrix.
**Symptom:** An ID-only matrix, or a suggestion to use non-existent `Intensity <sample>` columns.
**Fix:** Use the corrected reporter columns shown above, retain the channel-to-sample map, and route reporter correction/normalization to quantification.

### Zero treated as a measurement

**Trigger:** `np.log2` applied directly to a MaxQuant matrix still containing 0.
**Mechanism:** MaxQuant encodes "not quantified" as 0; log2(0) = -inf, which then propagates into means and tests.
**Symptom:** -inf values, NaN means, proteins silently dropped or skewed.
**Fix:** `replace(0, np.nan)` before any transform; then diagnose missingness.

### Bookkeeping rows survive

**Trigger:** Loading proteinGroups.txt without filtering Reverse / Potential contaminant / Only identified by site.
**Mechanism:** Decoys exist only for FDR estimation; contaminants are keratin/trypsin/BSA, not the sample; site-only groups have no unmodified-peptide quant evidence. `Only identified by site` exists only in proteinGroups.txt.
**Symptom:** Inflated protein counts; a "hit" that is a decoy or keratin.
**Fix:** Filter all three flag columns; cross-check with `REV__`/`CON__` ID prefixes when joining to peptide tables. Caveat: do not delete CON__ rows blindly if a contaminant (e.g. keratin) is the protein of interest.

### Razor / leading protein-ID ambiguity ignored

**Trigger:** Treating `Protein IDs` or `Gene names` as an atomic single value.
**Mechanism:** These are semicolon-delimited lists; the first entry is the leading (razor) protein for the group, and `Gene names` can be blank while protein IDs are present.
**Symptom:** Merges fail, NaN gene labels, ambiguous identity downstream.
**Fix:** Split on `;` and take the first entry; guard `Gene names` with `.notna()`. Group parsimony details -> protein-inference.

### Stale DIA-NN parsing

**Trigger:** Reading `report.tsv` on DIA-NN 2.0, or pivoting before q-filtering.
**Mechanism:** 2.0 defaults to (and only defaults to) `report.parquet`; pivoting unfiltered rows includes precursors above 1% FDR.
**Symptom:** FileNotFoundError on report.tsv; or low-confidence quant inflating the matrix.
**Fix:** `pd.read_parquet('report.parquet')`; filter `Q.Value <= 0.01 & PG.Q.Value <= 0.01 & Global.PG.Q.Value <= 0.01` before pivoting `PG.MaxLFQ`.

### MNAR imputed as MCAR

**Trigger:** Mean/median/KNN imputation on a DDA or DIA matrix whose diagnostic shows abundance-dependent missingness.
**Mechanism:** DDA missingness is abundance-dependent (left-censored); MCAR imputers fill missing low values with the central tendency, biasing them upward.
**Symptom:** Low-abundance proteins gain false high values; spurious differential hits.
**Fix:** Diagnose the abundance-missingness correlation here and choose the imputer from it, not from the acquisition mode: a negative correlation routes to left-censored imputation (downshifted-Gaussian / QRILC / MinProb) or a censoring-aware model in differential-abundance. DIA's fewer missing values limit the damage but do not make them MCAR.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| DIA-NN import filter `Q.Value <= 0.01` AND `PG.Q.Value <= 0.01` AND `Global.PG.Q.Value <= 0.01` | Demichev 2020; DIA-NN documentation | Run-level precursor/protein and experiment-wide protein-group 1% FDR enforced before any quant value enters a cross-run matrix |
| Peptide/protein FDR 1% (q <= 0.01) | Target-decoy convention | Standard ID confidence at both peptide and protein levels |
| MaxQuant zero -> NaN | MaxQuant output convention | 0 encodes "not quantified"; log2(0) = -inf corrupts every transform |
| Min peptides per protein for quant >= 2 | Community quant practice | Single-peptide ("one-hit-wonder") proteins are ID/quant-unreliable; not applied by the import code above -- apply it (e.g. `Razor + unique peptides >= 2`) as a documented study choice |
| Valid-value filter >= 50-70% per group | Modeling choice (document per study) | Caps imputation burden; the exact cutoff is a study decision, not a universal constant |
| Take FIRST semicolon entry as leading protein/gene | MaxQuant proteinGroups convention | The leading/razor protein is the group identifier; trailing entries are shared-peptide members |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `-inf` values after log2 | Zeros not converted to NaN (MaxQuant LFQ or DIA-NN PG.MaxLFQ) | `df.replace(0, np.nan)` before `np.log2` |
| `IndexError: list index out of range` at `getPrecursors()[0]` | All-ion (AIF/MSE) MS2 scan without a precursor | Check `if not spectrum.getPrecursors()` before indexing |
| `FileNotFoundError: report.tsv` (DIA-NN 2.0) | TSV no longer the default output | `pd.read_parquet('report.parquet')` |
| `KeyError: True` / bookkeeping fields absent | A non-MaxQuant table was passed to the proteinGroups reader; scalar `.get(..., '')` masks are unsafe when every flag is absent | Build a row-indexed mask; fail clearly if all three MaxQuant bookkeeping fields are absent |
| No `LFQ intensity` but `Reporter intensity corrected` exists | This is a MaxQuant TMT table, not label-free data | Import corrected reporter channels; do not use the single `Intensity` total; reporter processing -> quantification |
| Mixed-type / DtypeWarning on MaxQuant load | Wide TSV with mixed column types | `pd.read_csv(..., low_memory=False)` |
| NaN gene labels break a merge | `Gene names` is a semicolon list, sometimes blank | `.where(notna(), '').str.split(';').str[0]` |
| Ratios track loading not biology | Read `Intensity` (raw) instead of `LFQ intensity` | Use `LFQ intensity` for between-sample comparison |
| `get_peaks()` unpacking error | Expecting a 2D array | It returns a tuple `(mz, intensity)` of two numpy arrays |

## References

- Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. 2014. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. *Mol Cell Proteomics* 13(9):2513-2526.
- Demichev V, Messner CB, Vernardis SI, Lilley KS, Ralser M. 2020. DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. *Nat Methods* 17(1):41-44.
- Rost HL, Schmitt U, Aebersold R, Malmstrom L. 2014. pyOpenMS: a Python-based interface to the OpenMS mass-spectrometry algorithm library. *Proteomics* 14(1):74-77.
- Hediyeh-zadeh S, Webb AI, Davis MJ. 2023. MsImpute: estimation of missing peptide intensity data in label-free quantitative mass spectrometry. *Mol Cell Proteomics* 22(8):100558.
- Chambers MC, Maclean B, Burke R, et al. 2012. A cross-platform toolkit for mass spectrometry and proteomics. *Nat Biotechnol* 30(10):918-920.
- Hulstaert N, Shofstahl J, Sachsenberg T, et al. 2020. ThermoRawFileParser: modular, scalable, and cross-platform RAW file conversion. *J Proteome Res* 19(1):537-542.

## Related Skills

- peptide-identification - search raw spectra and convert vendor RAW to mzML
- quantification - compute MaxLFQ and TMT reporter-ion quantities from imported data
- protein-inference - resolve protein-group parsimony and razor assignment
- differential-abundance - normalize, impute (per the missingness diagnosis), and test
- proteomics-qc - assess run-level identification and quant quality
- dia-analysis - run DIA-NN to produce the report this skill imports
- expression-matrix/normalization - general intensity-matrix normalization patterns
- workflows/proteomics-pipeline - end-to-end pipeline that begins with this import step
