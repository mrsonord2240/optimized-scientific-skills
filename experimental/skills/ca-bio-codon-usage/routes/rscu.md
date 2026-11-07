# Relative synonymous codon usage

```bash
python scripts/rscu_analysis.py
```

This runs the source's hard-coded illustrative CDS. For a supplied CDS, adapt the example or call `calculate_rscu(seq, table_id=..., policy='strict')`; no file-based interface is supplied. It divides observed counts by the uniform expectation within each synonymous family, reporting 0 for a family with no observations. The example labels observed RSCU below 0.5 as rare, with 0.8 and 1.2 as descriptive display thresholds. Select the correct genetic code and confirm the input frame and orientation. See [metrics](../references/metrics-and-methods.md).

**Done when:** RSCU values, genetic code, input and validation details are retained with the family-level interpretation.
