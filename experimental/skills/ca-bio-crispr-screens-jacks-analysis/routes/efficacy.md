# Summarize guide efficacy
~~~bash
python scripts/efficacy_summary.py jacks_out_grna_JACKS_results.txt guidemap.txt --low 0.3 --out low_eff_by_gene.tsv
~~~
Inputs require sgrna, X1, X2 and guide-map sgRNA, Gene. X1 is efficacy posterior mean; X2 is a second moment, not standard deviation. Threshold 0.3 is the source convention; missing map IDs raise an error. Review context and effect/std, not efficacy in isolation.
The source example scripts/analyze_results.py additionally plots fixed-prefix results and analyzes the first cell line; adapt it to the intended data.
Done when per-gene summary and low-efficacy limitations are reviewed. See [efficacy notes](../references/efficacy-prior-and-diagnostics.md).
