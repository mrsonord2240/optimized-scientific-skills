# Quantitative Thresholds: Effector-Gene Prioritization

Moved out of `SKILL.md` (bio-causal-genomics-effector-gene-prioritization) because the full
threshold table is needed only when scoring or reporting a specific evidence stream --
consult this file when `SKILL.md`'s Multi-Evidence Integration Framework or Reconciliation
sections point you here.

| Quantity | Threshold | Source / Rationale |
|----------|-----------|---------------------|
| L2G score (high-confidence) | >= 0.5 | Open Targets default; gradient-boosted classifier calibrated against curated gold standards |
| L2G score (suggestive) | >= 0.2 | Open Targets exploratory threshold |
| MAGMA gene-wide p | < 2.5e-6 (Bonferroni 0.05 / 20k genes) | Standard genome-wide gene-level significance |
| coloc PP.H4 (triangulation) | >= 0.7 | Open Targets / common practice; >= 0.8 for stringent |
| PoPS score (high-confidence) | Top decile per locus | Weeks 2023 Nat Genet 55:1267; threshold is relative per-locus rank, NOT an absolute cutoff. Absolute PoPS score is scale-dependent on trait polygenicity, so an absolute "PoPS >= 0.5" rule is incorrect across traits |
| ABC enhancer-gene score | >= 0.02 (standard) or >= 0.04 (stringent) | Fulco 2019; cross-reference atac-seq/enhancer-gene-linking |
| ENCODE-rE2G probability | >= 0.5 (binarised) | Gschwind 2023 |
| cS2G aggregate score | >= 0.5 per SNP-gene allocation | Gazal 2022; heritability-calibrated aggregator |
| Distance to TSS (regulatory window) | <= 100 kb (default); <= 500 kb (liberal); <= 1 Mb (absolute) | Convention; Mountjoy 2021. Beyond 100 kb distance ceases to be a reliable single feature |
| MAGMA gene-window | 35 kb upstream + 10 kb downstream | FUMA default; balances regulatory capture vs gene-dense dilution |
| Fine-mapping PIP (causal variant) | > 0.5 (suggestive); > 0.9 (strong) | Convention (cross-reference causal-genomics/fine-mapping) |
| Multi-evidence concordance | >= 3 of 6 streams | Operational rule from Open Targets Genetics and Mountjoy 2021 |
| Single-cell eQTL panel size | >= 200 donors per cell type | Below this, per-cell-type eQTL discovery underpowered |
