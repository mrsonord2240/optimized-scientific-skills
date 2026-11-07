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
| Correct QFeatures summary followed by Windows exit `3221225477` / `0xC0000005` | A native Windows dependency can crash during process teardown; one reproduced case was `cli.dll` in `cli__kill_thread`, after the analysis had finished | Treat the nonzero exit as failure. Reinstall/update `cli` in a clean R library, or run the same script in WSL/Linux; confirm a zero exit code rather than relying on printed output |
