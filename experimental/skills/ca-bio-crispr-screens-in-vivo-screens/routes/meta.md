# Combine per-sample results
~~~bash
python scripts/per_animal_meta_analysis.py
~~~
Run from the skill directory after adapting the fixed animal_files mapping to the input paths in the script. The source example reads up to ten MAGeCK gene_summary files; absent files are silently skipped. It combines negative-selection p-values by Stouffer Z and BH, writing in_vivo_meta_all.tsv and in_vivo_meta_hits.tsv. The source threshold is meta-FDR <0.05 with nominal p <0.05 in at least half of all loaded samples; meta-significant genes below that threshold remain in the all table for review.
Inspect missing genes, number of samples per gene and the shared baseline: ordinary Stouffer assumes independence, which shared-control comparisons may violate. The code clips only the lower p-value boundary, so p=1 can produce infinite z; check values and calibrate the analysis before interpreting results.
Done when sample coverage, dependence, thresholds and flagged cases are reviewed.
