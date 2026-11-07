# ArchR and SnapATAC2

```bash
Rscript scripts/archr_workflow.R samples.csv archr_run /path/to/macs3
```

`samples.csv` has `sample` and `fragments` columns; paths point to indexed fragment files. The wrapper creates a new `archr_run` directory and uses hg38, source QC thresholds, doublet filtering, TileMatrix LSI and cluster peak calling. For Python use `python scripts/snapatac_workflow.py fragments.tsv.gz output.h5ad`; this source workflow also assumes hg38 and runs MACS3. Both need their respective dependencies and suitable data depth. Adapt genome and annotations for other assemblies, and inspect empty Arrow results. Details and gene-activity examples: `references/ecosystem-workflows.md`.

Done when the requested project, clustering and relevant peak matrix are produced and QC recorded.
