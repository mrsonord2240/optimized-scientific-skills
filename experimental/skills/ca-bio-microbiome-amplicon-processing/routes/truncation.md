# Choosing truncLen

No script; the numbers go into `--trunc-len` of `scripts/dada2_asv.R` (pass `--amplicon-len` so it checks the budget).

truncLen is a merge-overlap budget first and a quality cut second. Paired reads merge only if `truncLen_F + truncLen_R >= amplicon_length + 12` (DADA2 `minOverlap` default 12). Cut where median quality falls below about Q25-30 on `plotQualityProfile`, and for long amplicons let the budget win.

Ceiling: truncLen must not exceed the read length after primer removal, not the raw cycle count. `filterAndTrim` drops, never pads, a shorter read; too high a value drops every read with only a `No reads passed the filter` warning. A 2x250 run with 19 bp 515F and 20 bp 806R removed leaves about 231 and 230 bp. Confirm with `plotQualityProfile` or `nchar` on a few trimmed reads. The script prints the trimmed read lengths and stops if truncLen exceeds the longest.

| Amplicon | Run | Slack | Do |
|----------|-----|-------|----|
| V4 515F/806R, about 253 bp | 2x250 | large | Truncate to quality within the post-primer length, e.g. 220,200 |
| V3-V4 341F/805R, about 460 bp | 2x250 | about 28 bp | Barely truncate the reverse; `--max-ee 2,5` to keep low-quality reverse reads |

`maxEE` is computed on the truncated read, so it interacts with truncLen.

Done when the sum of truncLen clears amplicon length + 12 and each value is at or below the trimmed read length.
