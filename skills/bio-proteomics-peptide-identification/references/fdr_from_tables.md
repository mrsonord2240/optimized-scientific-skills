# q-values from a results table

Read when you hold a PSM table from any engine (concatenated search) or separate target and decoy tables and need q-values. `examples/separate_search_fdr.py` ships the separate-search code as a script.

### FDR from a Results Table (concatenated competition, made explicit)

**Goal:** Compute q-values from any engine's PSM table when the search was a single concatenated target-decoy search.

**Approach:** Keep the best hit per spectrum, rank by score, walk down accumulating target and decoy counts, FDR = (decoys + 1)/targets, then take the running minimum from the bottom to get monotone q-values. `score` must be higher-is-better, and the decoy prefix must match the engine's (Sage writes lowercase `rev_`); a table with no recognised decoys must stop, not pass every PSM. This form is correct ONLY for concatenated competition; separate searches need pi0 * decoys/targets (Kall et al. 2008; pi0 = 1 is the conservative default) or the mix-max estimator (Keich, Kertesz-Farkas & Noble 2015; Percolator's default for separate-search input).

```bash
python scripts/table_fdr.py search_results.tsv --scan scan --score score --protein protein --out psms_1pct.tsv
```

`scripts/table_fdr.py` takes the three column names, prefix-matches the decoys lower-cased (`decoy_`, `rev_`, `xxx_`), raises when none are recognised (the message states the smallest reachable q, 1/rows), keeps the best hit per spectrum, and applies (decoys + 1)/targets with the running minimum. `score` must be HIGHER-is-better: Sage `sage_discriminant_score` or Comet `xcorr` as is; E-values (Comet `e-value`, MS-GF+ `SpecEValue`) as -log10(E-value), or it keeps nothing.

### FDR from SEPARATE Target and Decoy Searches (pi0 * D / T)

**Goal:** Get a valid 1% list out of two result tables produced by searching the same spectra against a target DB and a decoy DB independently.

**Approach:** No competition resolved which hit wins, so the decoy count estimates the number of incorrect TARGETS directly, scaled by pi0, the proportion of target PSMs that are incorrect (Kall et al. 2008). pi0 = 1 is always valid and conservative; the median-decoy estimate (twice the fraction of target scores below the median decoy score) recovers the identifications that pi0 = 1 throws away, at the cost of estimating a nuisance parameter. Do NOT run the concatenated-search script above on the two tables merged. `examples/separate_search_fdr.py` ships this as a script.

```bash
python examples/separate_search_fdr.py target.tsv decoy.tsv     # pi0-hat list and the pi0 = 1 bound
```

The functions `estimate_pi0` (median-decoy estimator, Kall et al. 2008) and `separate_search_qvalues(targets, decoys, pi0=None)` can be imported from that script; it also rejects duplicated scans per table and an empty decoy table.

On the synthetic separate-search pair with ground truth (12,000 spectra each): pi0-hat = 0.611 keeps 2,888 PSMs at a true FDP of 1.04%, pi0 = 1 keeps 2,588 at 0.62%, and Elias-Gygi's 2d/(t+d) misapplied here keeps only 2,139 at 0.33%. The alternative is to hand both tables to Percolator and let mix-max do it: that is Percolator's default for separate-search input, and it is a calibrated-score procedure, not the same arithmetic.
