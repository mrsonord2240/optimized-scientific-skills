# TMT channel balance

```bash
python scripts/tmt_balance.py plexA.csv plexB.csv --out tmt_balance_results.tsv
```

One CSV per plex: protein id in the first column, one column per channel, RAW reporter intensities (before any channel or global normalization). Plex name is the file stem. The script refuses a table that looks log-transformed or normalized.

- Each channel is compared with its own plex median. Plexes differ in overall signal, so never compare channels across plexes; use the reference channel for that.
- `investigate` is beyond `--max-fold` (default 2x) from the plex median; `flag_severe` beyond 3.5x. Each investigated channel is a row in the exclusion log (`routes/writeup.md`).
- MSstatsTMT `dataProcessPlotsTMT` shows balance only after `proteinSummarization(..., global_norm = FALSE, reference_norm = FALSE)`; its default equalizes channel medians and hides imbalance. It needs PSM-level input.
- Wide isolation windows compress every fold change toward 1 (a real 10:1 reads about 5:1). Filter isolation interference below 50% (below 30% is stricter), prefer SPS-MS3, FAIMS or narrow windows, and measure the floor with an empty or TKO channel.

Done when every channel's fold against its plex median is in `tmt_balance_results.tsv` and the investigated ones are named.
