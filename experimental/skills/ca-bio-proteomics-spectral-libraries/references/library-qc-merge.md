# Library QC and merge

Read when summarising a library or merging several libraries.

## QC and Merge Libraries

**Goal:** Summarize a library and combine multiple libraries without dropping legitimate distinct transitions.

**Approach:** Report precursor/protein counts and transitions-per-precursor, then dedup on the FULL transition key. Deduping on (sequence, fragment-type, fragment-number) alone drops real transitions that differ only in precursor charge or fragment charge -- key on all five.

`merge_libraries` and `library_stats` (with `TRANSITION_KEY`) are exposed by the file-based merge route in `scripts/library_tools.py`. The key is `ModifiedSequence, PrecursorCharge, FragmentType, FragmentSeriesNumber, FragmentCharge`.
