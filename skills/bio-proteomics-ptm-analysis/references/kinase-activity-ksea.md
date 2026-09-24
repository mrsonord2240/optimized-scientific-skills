# Kinase Activity with KSEAapp

**Goal:** Score kinase activity from the protein-adjusted site fold-changes.

**Approach:** `KSEAapp::KSEA.Scores` merges the site table with a kinase-substrate prior on gene symbol + residue and returns a z-score per kinase with its substrate count `m`. The prior (`KSData`: PhosphoSitePlus + NetworKIN table with `KINASE`, `SUB_GENE`, `SUB_MOD_RSD`, `Source`, `networkin_score` columns) must be the full file from github.com/casecpb/KSEA; the package's `data(KSData)` is an abbreviated demonstration subset. `PX` needs six columns in this exact order, and `FC` is a LINEAR ratio (the function takes log2 itself), treatment over control.

**Input boundary:** Use only the checked runner's published paired-global label-free or TMT result for regulatory kinase-activity calls. In the checked examples this is the versioned path held in `FINAL`; any retained private staged `adjusted_sites.csv` / `adjusted_sites_tmt.csv` is evidence only and is not the published result. Do not pass `proxy_adjusted_sites.csv` or `proxy_ptm_model.csv` from a no-global `use_unmod=TRUE` analysis into KSEA: those are proxy-adjusted candidates, not protein-adjusted regulation evidence.

```bash
Rscript scripts/ksea_scores.R adjusted="$FINAL" proteinGroups=proteinGroups_global.txt prior='PSP&NetworKIN_Kinase_Substrate_Dataset.csv' out=kinase_scores.csv
```

`scripts/ksea_scores.R` builds PX from the adjusted table (`FC = 2^log2FC`, gene symbols from the global `proteinGroups.txt`), drops non-finite `log2FC` rows before PX, stops with a message when no finite site, no gene symbol or fewer than 2 prior-covered sites remain, then calls `KSEA.Scores(KSData, PX, NetworKIN = FALSE, NetworKIN.cutoff = 3)` and writes `Kinase.Gene, m, z.score, FDR` sorted by z-score.

If the contrast Label reads `Control vs Treatment`, use `FC = 2^(-log2FC)`, or every kinase's sign inverts.
