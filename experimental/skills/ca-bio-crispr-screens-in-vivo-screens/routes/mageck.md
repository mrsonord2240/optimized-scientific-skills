# Analyze existing count data
~~~bash
mageck test -k in_vivo_screen.count.txt -t Animal1 -c Plasmid -n animal_1
~~~
Run the same source-specified comparison separately for each sample, supplying exact count-table column names and distinct output prefixes. This produces per-sample gene_summary tables for [meta-analysis](meta.md). Preserve sample provenance and assess shared-control dependence. For FASTQ-to-count use bio-crispr-screens-mageck-analysis and check installed version.
The source also proposes an MLE design with one plasmid and one observation per animal, using animal indicators plus a tumor indicator. That example is saturated: it leaves no residual degrees of freedom and does not separate replicate variation from condition; do not treat it as a validated batch-adjusted model. A defensible model requires a suitable replicated design and independent review.
Done when count/sample mappings, model assumptions and per-sample output quality are reviewed.
