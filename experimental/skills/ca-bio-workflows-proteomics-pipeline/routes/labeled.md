# Labeled proteomics
~~~bash
Rscript scripts/tmt_impurity_correct.R tmt.mzML tmt10_coa.csv reporters.csv
Rscript scripts/msstatstmt_multiplex.R evidence.txt proteinGroups.txt annotation.csv multiplex.csv
Rscript scripts/silac_limma.R proteinGroups.txt silac.csv
~~~
Choose the relevant route and a new output directory; scripts may overwrite. Done when TMT impurity orientation, per-plex Norm bridge and the output are checked, or SILAC ratio quality and labeling are reviewed. The multiplex source checks only that Norm exists somewhere: independently verify it exists in every plex. TMT CoA rows must be source channels and columns observed channels. SILAC has Arg-to-Pro and labeling-efficiency caveats. See [TMT](../references/tmt-isobaric.md) and [SILAC](../references/silac.md).
