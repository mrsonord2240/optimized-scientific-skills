# Inspect and preserve artifacts

```bash
qiime tools peek table.qza
```

Use `qiime tools validate table.qza --level max` to check archive and format. `qiime tools extract --input-path table.qza --output-path extracted` preserves the archive structure for inspection; `qiime tools export --input-path table.qza --output-path exported` extracts native data without the artifact's provenance. Retain the original artifact and record the export boundary. Replay commands and citations with `qiime tools replay-provenance --in-fp table.qza --out-fp replay.sh --usage-driver cli` and `qiime tools replay-citations --in-fp table.qza --out-fp citations.bib` after checking installed help. Replay is not a promise of bit-identical results across releases. A visualization `.qzv` is terminal; retain its input `.qza`.

Done when type, validation, provenance availability, release and any export boundary are documented.
