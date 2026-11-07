# Python only

```bash
python scripts/welch_de.py proteinGroups.txt samples.csv Treatment-Control results.csv
```

- Per-protein Welch t-test and Benjamini-Hochberg after log2 and median centring. No variance moderation.
- The script warns below 10 samples per group. At 3 to 5 per group the p-values are unreliable: say the result is exploratory, or install R and use `routes/protein-table.md`.
- Other column prefix: `prefix="LFQ intensity "`. A CSV of linear intensities (first column protein ID) also works.
- Proteins with fewer than `min_obs=2` values in a group go to `results_untestable.csv`. Report them as undetected.
- There is no Python equivalent of ashr; shrunken fold changes need R (`routes/fold-change.md`).

Done when `results.csv` and `results_untestable.csv` exist.
