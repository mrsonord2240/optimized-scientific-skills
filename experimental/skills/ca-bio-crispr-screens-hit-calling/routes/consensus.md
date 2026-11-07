# Compare method calls
~~~bash
python scripts/consensus_hits.py mageck.gene_summary.txt bagel.bf.txt drugz.txt -o consensus.tsv
~~~
The three inputs must represent the same biological contrast and compatible directions. MAGeCK needs id and neg|fdr, BAGEL2 GENE and BF, drugZ GENE and fdr_synth. Defaults are FDR <0.05 and BF >6. The script outer-joins genes and flags overlap not enriched above chance; a warning does not prove mismatched comparisons, and no warning does not establish comparability. Missing-method rows count as not called. For two-method MAGeCK/BAGEL2, run scripts/consensus_two.py with fixed mageck.gene_summary.txt and bagel_bf.txt; it writes consensus_hits.csv.
Done when comparison identity, missingness, method assumptions and consensus tiers are reviewed. See [failure modes](../references/failure-modes.md).
