# Evidence framework
When fine-mapping returns credible sets (cross-reference causal-genomics/fine-mapping), each credible variant should be assigned to a gene using a fixed lexicographic ladder rather than a single feature. The ladder collapses ambiguity by preferring direct mechanistic evidence first and falling back to weaker signals only when stronger ones are absent:

1. **Coding consequence at credible variant** (missense, splice-donor / splice-acceptor, stop-gained, start-lost via VEP / Ensembl consequence) -> assign variant to that gene.
2. **eQTL / pQTL colocalization PP.H4 >= 0.7** with the gene's expression QTL (matched tissue) -> assign to that gene (cross-reference causal-genomics/colocalization-analysis).
3. **ABC or ENCODE-rE2G enhancer-gene linkage** when matched ATAC + H3K27ac (+ optional Hi-C) is available in the candidate tissue -> assign to the linked gene (cross-reference atac-seq/enhancer-gene-linking).
4. **Nearest TSS** -> assign as last-resort fallback; flag as low-confidence (nearest-gene is correct only 50-70% of the time at well-fine-mapped loci).

If a single credible variant ties across two or more genes at the same rung (e.g. coding consequence in gene A AND a competing eQTL coloc to gene B), report multi-gene candidacy explicitly rather than forcing a single assignment; the locus may be multi-effector or the credible set may straddle a regulatory boundary.

## Multi-Evidence Integration Framework

A strong candidate causal gene at a GWAS locus requires concordance across multiple orthogonal evidence streams. The six canonical streams:

| Evidence stream | What it tests | Pass threshold | Source |
|----------------|---------------|----------------|--------|
| Fine-mapping | SuSiE PIP for variant; variant assigned to gene by ABC or eQTL coloc | PIP > 0.5; credible set purity > 0.5 | susieR (cross-reference causal-genomics/fine-mapping) |
| Colocalization | Shared causal variant with gene's eQTL / pQTL / sQTL | coloc.abf or coloc.susie PP.H4 >= 0.7 | coloc (cross-reference causal-genomics/colocalization-analysis) |
| Distance | Variant within annotated regulatory unit (gene body or enhancer-gene unit) | Distance to TSS <= 100 kb OR within ABC enhancer-gene unit | Convention; Mountjoy 2021 |
| Polygenic prior (similarity) | Gene shares pathway / co-expression / PPI with other GWAS hits for the trait | PoPS score in top decile per locus | Weeks 2023 |
| L2G (integrative classifier) | Per-locus per-gene gradient-boosting on a panel of features | L2G score >= 0.5 (Open Targets default high-confidence) | Mountjoy 2021 |
| Chromatin / enhancer-gene | ABC or ENCODE-rE2G connects fine-mapped variant to gene's promoter | ABC >= 0.02 OR ENCODE-rE2G >= 0.5 | Fulco 2019; Gschwind 2023 |
| Deep-learning variant effect (optional 7th) | chromBPNet / EnFormer in silico variant effect on accessibility/expression at credible variant | \|log2FC\| > 1 (chromBPNet strong-effect); agreement across two models | chromBPNet: Pampari 2025 bioRxiv 2024.12.25.630221 (preprint); EnFormer: Avsec 2021 Nat Methods 18:1196; cross-reference atac-seq/deep-learning-atac |
| Cross-trait coincidence (optional 8th) | Same gene flagged at related-trait loci (e.g. lipid GWAS at CHD lead variant) | Same gene top-ranked at >= 2 related traits | Convention |

**Operational rule:** Report a gene as a "high-confidence causal effector" only when >= 3 of the 6 core evidence streams are concordant (deep-learning variant effect and cross-trait coincidence are optional supplementary streams). >= 4 concordant is "strong-confidence"; >= 5 is "near-certain pending experimental validation". Single-stream evidence is associational only; two-stream concordance is suggestive. This is the standard used by Open Targets Genetics (Mountjoy 2021), GTEx-derived target nomination pipelines (Open Targets Platform), and pharma drug-discovery workflows.

## Quantitative Thresholds

Canonical pass/fail thresholds for every evidence stream (L2G, MAGMA gene-wide p, coloc PP.H4,
PoPS decile, ABC/ENCODE-rE2G, cS2G, distance-to-TSS, MAGMA window, fine-mapping PIP,
multi-evidence concordance, single-cell eQTL panel size) live in
`references/quantitative-thresholds.md`. Consult it when scoring or reporting a specific
evidence stream -- the two thresholds used constantly are repeated inline where needed
(coloc PP.H4 >= 0.7 in the Multi-Evidence table below; >= 3 of 6 concordance in the
Operational rule).
