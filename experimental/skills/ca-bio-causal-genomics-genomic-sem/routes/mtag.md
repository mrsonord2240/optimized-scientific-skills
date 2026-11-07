# MTAG comparison

```bash
bash scripts/mtag_pipeline.sh
```

Edit trait input paths first and set MTAG_BIN and output location. Inspect per-trait MaxFDR; the source treats > 0.05 as unreliable. MTAG produces per-trait effects rather than common-factor effects; compare as described in `references/mtag-comparison.md`.

Done when trait-specific results, MaxFDR and any discordance are reported.
