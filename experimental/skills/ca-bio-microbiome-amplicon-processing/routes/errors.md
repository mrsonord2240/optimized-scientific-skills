# Errors and wrong-looking results

| Symptom | Cause | Fix |
|---------|-------|-----|
| Near-zero merge rate | truncLen below the overlap budget | Recompute the budget; keep length, loosen reverse `maxEE` (`routes/truncation.md`) |
| Large read fraction called chimeric | Primers not trimmed (degenerate bases) | `scripts/remove_primers.sh` before filtering |
| `plotErrors` fitted line diverges from points | Binned quality (NovaSeq, NextSeq) | Enforce a monotonic error matrix; nf-core/ampliseq `--illumina_novaseq` |
| ASVs vanish or appear when runs are split | One error model fit across runs | Per-run `learnErrors`, then `mergeSequenceTables` |
| Zero or few reads pass the filter | `maxEE` too strict, or truncLen above the post-primer read length | Loosen `maxEE`; set truncLen below the trimmed read length |
| ITS taxa lost, poor merging | Fixed `truncLen` on ITS | ITSxpress, then `--single --trunc-len 0` |

Done when the symptom's cause is fixed and the run repeated.
