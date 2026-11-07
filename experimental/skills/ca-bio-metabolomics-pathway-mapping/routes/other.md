# Raw features and network diffusion

```bash
Rscript scripts/fella_diffusion.R cpd_ids.txt fella_hsa fella_results.csv
```

Use KEGG compound IDs, one per line. On first use this source script builds a human KEGG graph via the network and caches it; inspect excluded compounds, resource use and output. Its default is normality-approximation diffusion, threshold 0.05. Mummichog/PSEA requires the entire feature table, correct ionization mode and mass tolerance. The retained `scripts/mummichog_psea.R` sets `mpt` input without RT while choosing v2, which the source says requires RT; it is not a ready-to-run path until the input format and version-specific API are reconciled. Ranked MSEA/QEA likewise needs input-specific adaptation; its MetaboAnalystR API route may transmit user data. Disclose external calls before use.

Done when method, full background, mapping exclusions, database and uncertainty are documented; network links and enrichment are hypotheses, not measured flux.
