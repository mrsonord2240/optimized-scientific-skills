# Compare continuous-trait models
```bash
Rscript scripts/continuous_file.R species_tree.nwk traits.csv > continuous-results.txt
```
**Done when:** model comparison and supported ancestral results with uncertainty are captured, or the script explicitly stops for an unsupported reconstruction. CSV row names must match tree tips exactly and include numeric `mass`; edit the data schema only after review. The source compares BM/OU/EB/lambda/kappa/delta. BM uses fastAnc confidence intervals; OUwie results are exploratory point estimates without uncertainty; other winning models stop rather than substituting BM. Review polytomy resolution and phylogenetic signal before interpretation. See [method workflows](../references/method-workflows.md).
