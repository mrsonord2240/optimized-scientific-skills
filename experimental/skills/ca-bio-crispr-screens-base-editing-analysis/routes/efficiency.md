# Editing efficiency

```bash
python scripts/filter_by_editing_efficiency.py crispresso_outputs --target-pos 5 --target-base C --threshold 0.5 --out efficiency.tsv
```

Position is 1-indexed within the quantification window, not the amplicon. Set the base and threshold for the actual assay. Input fractions use the source-checked CRISPResso2 table schema; missing sample tables are skipped. See `references/screen-analysis.md`, `references/crispresso-examples.md` and `references/be-validation-pipeline.md`.

Done when efficiency, missing samples and assay-specific filtering decisions are recorded.
