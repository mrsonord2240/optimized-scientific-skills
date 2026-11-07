# CEGv2 Essentialome Recovery

## Essentialome Recovery (CEGv2 PR-AUC)

**Goal:** Verify the screen has detectable biological essentiality signal by checking whether known essentials (Hart 2017 CEGv2) drop out faster than known non-essentials (NEGv1).

**Approach:** Compute precision-recall AUC where positives are CEGv2 genes and negatives are NEGv1; the screen "passes" if PR-AUC >0.7 (community convention; the CEGv2/NEGv1 sets come from Hart 2017 / Hart 2014).

```bash
python scripts/essentialome_recovery.py gene_lfc.tsv CEGv2.txt NEGv1.txt   # gene_lfc.tsv: columns gene,lfc
```

`scripts/essentialome_recovery.py` (also importable: `essentialome_recovery(gene_lfc_df, cegv2_set, negv1_set)`) keeps only CEGv2 and NEGv1 genes, scores by -LFC (negative LFC = depleted = more essential) and reports `pr_auc`, `roc_auc` and the essential / non-essential gene counts detected.

**Source / threshold:** Hart 2017 *G3* 7:2719 defines CEGv2 (~684 core essentials); NEGv1 (~927 non-essentials) comes from Hart 2014 *Mol Syst Biol* 10:733. Both lists at https://github.com/hart-lab/bagel/blob/master/CEGv2.txt and NEGv1.txt. DepMap convention: PR-AUC >0.7 at FDR 5% is the "passing" threshold; <0.5 means the screen has no essentiality signal and is not interpretable.

**When PR-AUC is low despite good Gini and Pearson**: cause is usually one of (a) Cas9 was not selected for before screen start (lots of Cas9-negative cells in the pool diluting signal), (b) puromycin selection truncated too aggressively (over-bottleneck), (c) the timepoint is too early (need 14-21 days for KO + decay + selection to manifest). Each has a different remediation.
