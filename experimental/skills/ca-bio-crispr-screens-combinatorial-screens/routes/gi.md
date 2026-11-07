# Score paired interactions
~~~bash
python scripts/gi_scoring.py
~~~
Run from a working directory with paired_lfc.tsv (cassette_id, gene_A, gene_B, lfc) and single_lfc.tsv (gene, lfc); adjust the script path when outside the skill directory. It writes gi_scores.tsv, synthetic_lethal_pairs.tsv and synthetic_rescue_pairs.tsv. It drops missing singleton pairs. GI equals observed double LFC minus the singleton sum. For small sets rank raw GI; z-scores have unstable nulls below roughly 50 pairs. Genome-scale calls use two-sided normal-tail p-values and BH 0.05, which rely on their null assumptions. Input duplicate singleton genes overwrite earlier values; check uniqueness and cell-line matching. Inspect output and validate independently.
Done when outputs, omissions and null assumptions are reviewed. See [scoring guidance](../references/scoring.md).
