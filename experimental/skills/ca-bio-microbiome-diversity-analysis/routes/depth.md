# Select sampling depth

```bash
Rscript scripts/pick_sampling_depth.R counts.tsv 0.01 0.85 40
```

Input is a feature-by-sample integer count TSV with feature IDs in the first column. The source script computes expected rarefied richness over a grid, seeks a relative gain below 1% per 1,000 reads among depths retaining at least 85% of samples, prints dropped samples and writes `depth_curve.tsv` in the working directory. If no plateau is found it reports the deepest allowed depth instead; this is a heuristic, not proof of saturation. Check counts and grid assumptions, assess excluded sample characteristics, and compare nearby depths. Output may be overwritten.

Done when depth, retained and dropped IDs, curve behavior and sensitivity are documented.
