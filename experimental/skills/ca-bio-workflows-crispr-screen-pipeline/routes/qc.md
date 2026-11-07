# QC of a guide count table

```bash
python scripts/qc.py experiment.count.txt qc.tsv plasmid=Plasmid
```

Column 1 is the guide id, column 2 the gene, then one raw integer column per sample. Replicates are grouped by stripping a trailing `_r1`, `_rep2`, `_A`, `_2` or a letter after a digit from the sample name; otherwise pass `pattern=<regex>`. `min_depth=300` is the reads-per-guide gate.

Without `plasmid=` the script warns and holds every sample to the endpoint Gini gate. The script prints a per-sample table and the within-condition replicate Pearson, writes `qc.tsv`, and exits 1 on a failed gate.

Gates:

- Gini on ln(count+1), the definition `mageck count` reports: plasmid below 0.1, endpoint below 0.2. `plasmid=` names the cloned library pool, not Day-0 cells.
- Reported, not gated: the share of guides above 25 reads (aim for 99% in the plasmid) and skew (p90/p10). The staged real plasmids give 99.0% and 4.3 (HAP1 TKOv3), 96.1% and 7.9 (A375 KY library), so judge them against the library, not a fixed cutoff.
- Replicate Pearson on log-counts at least 0.8, the MAGeCK-VISPR guideline. If no sample names group into replicates the verdict is `QC INCOMPLETE` (exit 1): rerun with `pattern=`.
- At least 300 reads per guide per sample.
- CEGv2 PR-AUC above 0.7 once endpoint-vs-baseline fold changes exist: `BAGEL.py fc` then `pr` (`routes/bagel2.md`), or the screen-qc Skill. Below 0.7 the screen failed selection; between 0.5 and 0.7, tighten FDR to 0.01.
- Average only within-condition replicate pairs: averaging every pair mixes baseline against endpoint and understates concordance.
- A failing replicate: compare each replicate's pairwise Pearson; if one replicate is clearly lower than the others, it is the outlier, so state that and rerun without it. If all pairs are similarly low, there is no outlier: report the screen's replicate concordance as failed.

Done when every gate is reported. A failed gate is stated in the answer, not hidden.
