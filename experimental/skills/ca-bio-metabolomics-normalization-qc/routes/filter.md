# QC feature filter

```bash
Rscript scripts/robust_dratio_filter.R peaks.csv filtered.csv sample_type QC 0.5 0.3
```

CSV has samples as rows, numeric feature columns and one sample-type column. Remove other metadata columns before this script; non-QC rows are treated as biological samples, so exclude blanks and conditioning samples separately. It retains features with MAD-based QC/biological D-ratio <= 0.5 and QC RSD <= 0.3; undefined values are removed. These defaults are conventions, not universal requirements. Record stage, input counts, exclusions and thresholds; this combined filter does not separately report each criterion's removed count. Output may be overwritten.

Done when QC and biological sample labels, filtered feature counts and threshold suitability are checked.
