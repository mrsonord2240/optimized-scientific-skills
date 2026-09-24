<!-- Moved from SKILL.md (2026-09-21); code now in scripts/. -->

## Mummichog / PSEA on a Raw m/z Peak Table

**Goal:** Predict perturbed pathway activity from an untargeted LC-MS feature table when no compound identities exist.

**Approach:** Declare instrument ppm and ionization mode, load the FULL feature table (m/z + p-value + t-score, optionally RT), set the query-defining p-cutoff, and run PSEA whose permutation null is sampled from R_all.

```bash
# peaks.csv: the ENTIRE feature table (m/z, p-value, t-score), not just significant peaks
Rscript scripts/mummichog_psea.R peaks.csv 5.0 negative 0.2 1000 mummichog_psea.csv   # ppm, ionization mode, query p-cutoff, permNum
```

The script sets `set.seed(123)` (PerformPSEA's permNum resampling is not reproducible otherwise), uses format `mpt` (`mprt` adds RT; use with `v2`), and `SetPeakEnrichMethod(mSet, 'mum', 'v2')` (`mum`|`gsea`|`integ`). ppm and ionization mode are chemistry-specific and mandatory: pos and neg use different adduct tables, and mixed data needs a per-feature mode column. The peak table must be the ENTIRE feature table: the permutation null draws random feature lists from it (R_all), so significant-only input pre-enriches the pool and makes everything significant. The query p-cutoff default is NOT 0.05; document the value used. The library string (`hsa_mfn`) encodes organism+network. The result `mSet$mummi.resmat` is predicted-active pathways, NOT a metabolite ID list.

`Read.PeakListData`/`SanityCheckMummichogData` log lines such as `A total of 11 of duplicates were
merged` (four passes on the audit's 1500-feature table): duplicate m/z-matched features are merged
before enrichment, so the feature count PSEA reports can be lower than the input row count.
