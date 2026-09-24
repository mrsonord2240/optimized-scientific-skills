# Specialised Methods and Reconciling Disagreeing Methods

## When Standard Pipeline is Insufficient

The FUSION / S-PrediXcan / S-MultiXcan + FOCUS pipeline is the default. Escalate to a specialised method only when the standard pipeline misses an expected hit (a strong GWAS signal with no TWAS gene, or a known causal gene without recovery).

| Method | Triggering scenario | Yield |
|--------|--------------------|-------|
| MOSTWAS (Bhattacharya 2021 PLoS Genet 17:e1009398) | Trans-mediator architecture suspected (immune, neuropsych traits) | ~15% additional hits via distal mediator terms |
| EpiXcan (Zhang 2019 Nat Commun 10:3834) | Paired epigenome data available (Roadmap, EpiMap, ENCODE cell-matched DNase/H3K27ac) | Higher prediction R^2 in epigenome-rich tissues |
| TIGAR-V2 (Parrish 2022 HGG Adv 3:100068) | Training a custom Bayesian DPR panel (no pre-trained PredictDB for the tissue) | Captures non-elastic-net effect structures |
| kTWAS (Cao 2021 Brief Bioinform 22:bbaa270) | Rare-variant or population-specific contexts; cis-eQTL panels under-powered | Kernel aggregation robust to non-linear and rare effects |

**Operational rule:** Default to FUSION / S-PrediXcan / S-MultiXcan + FOCUS first; document an expected-but-missing hit before escalating.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| FUSION significant + S-PrediXcan null | Different weight panels (FUSION vs PredictDB) or different LD references | Re-run with matched panel; check FUSION weight CV-R^2 vs PredictDB MASHR posterior mean |
| S-MultiXcan p << min per-tissue p | Joint test borrowing strength across tissues | Genuine; report the joint p and the top contributing tissue |
| Many co-significant genes at one locus, no FOCUS PIP > 0.5 | LD-tied co-regulated genes, true causal gene not in panel | Expand panel (add brain or tissue-specific weights); functional follow-up needed |
| FOCUS PIP > 0.8 + coloc PP.H4 < 0.5 | Sparse eQTL signal (single SNP drives prediction) + locus has another co-localising signal | Investigate the single top-eQTL SNP; check for fine-mapped credible-set overlap |
| TWAS hit + cis-eQTL MR null | TWAS hit is LD-tagged, not mediated by expression | Trust the cis-eQTL MR result; flag the gene as TWAS-positive but non-causal |
| MA-FOCUS PIP much lower than per-ancestry FOCUS | Cross-ancestry heterogeneity; gene effect is not shared | Report per-ancestry separately; do not pool |
| TWAS significant in wrong tissue | LD-induced via tissue-shared eQTL | Verify with LDSC-SEG tissue prioritisation; treat top-tissue TWAS hit as the trustworthy one |

**Operational rule:** No single-method TWAS result is sufficient evidence of causal gene identity. The minimum reporting standard is (a) TWAS significance threshold met with multiple-testing correction; (b) FOCUS PIP >= 0.8 OR coloc PP.H4 >= 0.7; (c) explicit acknowledgement of tissue choice and ancestry of prediction weights. Triangulation with cis-eQTL MR or independent CRISPR/MPRA validation lifts a finding from "candidate" to "supported".
