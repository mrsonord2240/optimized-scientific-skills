# Correct an isotopologue vector

```bash
python scripts/correct_csv.py intensities.csv FORMULA 0.99 corrected.csv
```

CSV must have ordered `intensity` values for M+0 through M+n; replace formula and measured per-position labeled purity. The wrapper invokes the source IsoCor low-resolution 13C helper and writes corrected areas and fractions. It does not infer formula, atom count, derivatization, high-resolution settings or data quality; verify those against the experiment and installed API. Check controls and corrected vector before interpretation. The helper discards IsoCor residuals; obtain them separately if needed for fit review. The source example has an approximate carbon-only numpy fallback; it is not a general replacement for full formula-aware correction and is not used by this wrapper. Source AccuCor high-resolution and multi-tracer methods require input-specific adaptation; set an owned output path. Output may be overwritten.

Done when formula, measured purity, method, controls and corrected MID are documented and checked.
