# Compare PE and BE screen results
~~~bash
python scripts/crossvalidate_pe_be.py be_screen_hits.tsv pe_screen_hits.tsv --fdr 0.05 --out pe_be_concordance.tsv
~~~
Input TSVs require variant_id and be_fdr/be_lfc or pe_fdr/pe_lfc. The script inner-joins variants and flags both FDRs below 0.05 with matching signs; unmatched variants are omitted and duplicate IDs can multiply rows. Inspect missingness, zero effects, comparability and edit coverage before interpretation.
Done when overlap and excluded variants are reviewed.
