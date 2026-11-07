# Filter prediction summaries
~~~bash
python scripts/filter_pridict2_summary.py predictions/summary.csv --threshold 50 --top 3 --out peg_library_filtered.csv
~~~
Use an existing PRIDICT2 batch summary with sequence_name and PRIDICT2_0_editing_Score_deep_K562 on a 0-100 scale; change --column for the corresponding HEK output. The >50 cutoff is a source project choice, not a prescribed model threshold. Report variants missing from results and prediction provenance. Source PRIDICT2 batch expects sequence_name and editseq, an existing fresh output directory, --cores 3, and --summarize with a supported model value; an empty summary may still exit zero. See [batch interface](../references/batch.md).
For existing edit-sequence input in the external PRIDICT2 environment, the source batch interface is:
~~~bash
python pridict2_pegRNA_design.py batch --input-fname variants.csv --output-dir predictions --cores 3 --summarize "$MODEL"
~~~
Set MODEL to a supported value and validate local help and input/output directories before use.
Done when column identity, units, missing designs and output are reviewed.
