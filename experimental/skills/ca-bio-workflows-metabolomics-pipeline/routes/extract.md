# Extract or import features
~~~bash
Rscript scripts/extract_cli.R metadata.csv positive extracted.rds
Rscript scripts/msdial_cli.R alignment.mdalign imported.rds
~~~
Choose one front end; run from the Skill root. Done when the new RDS contains a features-by-samples matrix and aligned metadata and extracted peaks have been inspected. CSV requires sample paths, sample_group, injection_order and batch in matching order. Acquisition mode is explicit; run modes separately. Filled peaks are imputations. MS-DIAL import excludes blank/standard rows; blank filtering requires a separate step. See [xcms](../references/stage1-xcms-extraction.md) and [MS-DIAL](../references/msdial-front-end.md).
