## MAJIQ V3 Workflow

**Goal:** Quantify LSVs with Bayesian posterior PSI distributions; ideal for complex multi-junction events that don't fit canonical event types.

**Approach:** Build a splice graph from BAMs + GFF3, compute per-junction coverage with bootstrap, then run `majiq psi` for posterior PSI per LSV.

**Not executed here:** MAJIQ V3 and VOILA need a registered academic/commercial licence (majiq.biociphers.org); the commands below and the V3 notes are unverified against a live install.

```bash
majiq build annotation.gff3 -c settings.ini -j 8 -o build_output
majiq psi build_output/sample1.majiq build_output/sample2.majiq -j 4 -o psi_output -n condition_psi
voila view -p 5000 -j 8 build_output/splicegraph.zarr psi_output/condition_psi.psi.voila -o voila_output
```

MAJIQ V3 (Aicher, Slaff, Jewell, Barash *bioRxiv* 2024; public release 2025) replaced V2's SQLite splicegraph (`splicegraph.sql`) with **Zarr storage** (`splicegraph.zarr`); the `.sql` is deprecated. V3 is substantially faster than V2 via a rewritten xarray/zarr implementation with parallelized coverage calculations. LSV output includes posterior mean PSI plus the full posterior distribution; this enables threshold-based testing (e.g. P(|ΔPSI| > 0.2)) rather than frequentist p-values.
