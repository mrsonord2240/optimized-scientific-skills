# Label-free MaxQuant analysis
~~~bash
Rscript scripts/limma_pipeline.R proteinGroups.txt sample_annotation.csv results.csv 0.6 1.5
~~~
Done when raw load QC, group completeness and design are reviewed and the protein-by-contrast output is inspected. Use a new output directory: source scripts can overwrite outputs. Annotation matches intensity suffixes and includes sample, condition, replicate and batch; first factor level is reference. The script drops contaminant flags indiscriminately: review target contaminants before use. It uses observed values and treat, not an explicit MNAR model; report on/off proteins as undetected. Median normalization assumes mostly symmetric changes; batch is a model covariate. See [guidance](../references/guidance.md).
