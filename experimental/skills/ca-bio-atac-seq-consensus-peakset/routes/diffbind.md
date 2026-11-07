# DiffBind summit re-centering

```r
library(DiffBind)
dba_obj <- dba(sampleSheet = 'samples.csv')
dba_obj <- dba.count(dba_obj, summits = 250, minOverlap = 2)
```

`summits = 250` re-centers every peak on its summit and gives 501 bp, matching the Corces width. The `dba.count` default is `summits = 200` (401 bp), so set it explicitly. This code is shown, not scripted; it needs MACS-style narrowPeak summit information.

- `minOverlap` is a count of samples. With unequal replicate numbers per condition a global majority rule drops a peak present in every replicate of the smaller condition: use `routes/per-condition.md` instead.
- Broad peaks lose width information under summit re-centering.
- DiffBind counting is N samples by M peaks and slow; set `DBA$config$cores` and `DBA$config$RunParallel`.
- DiffBind's consensus is denser than iterative overlap (summits within 500 bp are kept); both are valid.

Done when the count object exists and the answer states `summits` and `minOverlap`.
