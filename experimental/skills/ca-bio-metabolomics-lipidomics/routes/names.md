# Names and import

```bash
Rscript scripts/lipidr_import_checks.R lipids.csv converted.csv
```

CSV requires `Molecule` and sample intensity columns. The source wrapper converts sphingoid `;O1/2/3` suffixes to older lipidr prefixes, imports and warns about unclassified entries; review conversion and output before analysis. `python scripts/honest_level.py "LIPID_NAME"` parses supplied names and limits reported detail to molecular-species level at most; replace the placeholder and check each result against evidence. Goslin parsing does not validate the analytical assignment. Output may overwrite an existing file. See [interpretation](../references/interpretation.md).

Done when names, class failures, nomenclature and supported structural resolution are recorded.
