# Check what filtering did

```bash
bash scripts/filter_stats.sh before.vcf.gz after.vcf.gz wgs
```

Prints record counts, Ti/Tv before and after, and counts per FILTER label. The third argument (`wgs` or `wes`) sets the Ti/Tv floor for the warning. It warns when nothing was removed, when Ti/Tv fell after filtering (the filters are removing true transitions: relax them), and when Ti/Tv is below the floor (too loose).

| Metric | WGS | WES | Reading |
|--------|-----|-----|---------|
| Ti/Tv | 2.0-2.1 | 3.0-3.3 | Below range means excess false positives (random errors sit near 0.5). A WES set at about 2.1 is too loose |
| Het/Hom | 1.5-2.0 | 1.5-2.0 | Depends on ancestry. High suggests contamination, low suggests inbreeding or ROH. Stratify by ancestry before flagging |
| Known (dbSNP) fraction | above 99% | above 99% | Low recovery means over-filtering |

- A filter that improves one metric while degrading another is miscalibrated.
- Benchmark against a truth set with `bcftools isec`, stratified by GIAB BEDs, for an honest precision/recall.
- Full QC panel (het/hom by ancestry, contamination, relatedness): vcf-statistics Skill.

Done when before/after counts and Ti/Tv are reported with the assay named.
