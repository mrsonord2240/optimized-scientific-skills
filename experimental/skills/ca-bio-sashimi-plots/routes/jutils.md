# Cross-tool visualization

```bash
bash scripts/jutils_pipeline.sh /path/to/Jutils rmats_output/ meta.tsv bam_list.tsv annotation.gtf chr1:1000-2000 jutils_out/
```

This wrapper handles rMATS; converters for other inputs and metadata format are in `references/jutils.md`. Heatmap needs at least two passing events.

Done when the requested outputs are written and inspected.
