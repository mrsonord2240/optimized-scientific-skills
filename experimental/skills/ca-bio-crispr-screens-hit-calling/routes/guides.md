# Check guide support
~~~bash
python scripts/second_best_lfc.py mageck.sgrna_summary.txt --direction neg -o second_best.tsv
~~~
Input requires Gene and LFC columns. Inspect second-most-extreme LFC against a preselected threshold; single-guide genes return NaN and single_guide=True and cannot pass a second-guide check. Use --direction pos for enrichment.
Done when guide support and single-guide cases are reviewed; see [guide notes](../references/second-best-and-custom-zscore.md).
