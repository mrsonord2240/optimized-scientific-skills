# Compare doses
~~~bash
python scripts/dose_consistent_hits.py --top-dose high --out hits.tsv low=drugz_low.txt mid=drugz_mid.txt high=drugz_high.txt
~~~
Run drugZ separately for each dose against matched vehicle. Pass dose files in low-to-high order. The source rule requires a common normZ sign and top-dose directional FDR <0.05; monotonic effect size is supporting evidence, not required. Use --direction supp for suppressors; verify sign and direction yourself because the helper tests same sign separately from directional FDR and does not force sign to match direction. Missing genes are dropped; report them and inconsistent genes.
Done when consistent and omitted genes are reviewed; multi-factor designs may require bio-crispr-screens-mageck-analysis.
