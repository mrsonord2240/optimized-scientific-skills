# Variant score aggregation

```bash
python scripts/aggregate_variant_scores.py sgrna_summary.txt variant_annotation.tsv --out-prefix variant_scores
```

The annotation requires sgrna, target_variant and n_bystanders columns. The script separates target-only and mixed-guide outputs; reconcile editing efficiency and source sample identity before interpretation. See `references/screen-analysis.md` and `references/interpretation.md`.

Done when target-only and mixed outputs are recorded with evidence limits.
