# Method Comparison and Reconciliation

## Algorithmic Taxonomy

| Tool | Model | Input | Strength | Fails when |
|------|-------|-------|----------|------------|
| SuSiE / susie_rss (Wang 2020 JRSSB 82:1273; Zou 2022 PLoS Genet) | Iterative Bayesian sum-of-single-effects (IBSS), variational | Individual-level (X, y) or (z, R, n) | Fast; native PIP + credible sets; pluggable priors; default in modern pipelines | Reference LD mismatched to GWAS sample; locus dominated by polygenic background; >L true effects |
| SuSiE-inf / FINEMAP-inf (Cui 2024 Nat Genet 56:162) | SuSiE + infinitesimal random-effect component | (z, R, n) | Calibrated credible sets when locus is non-sparse (polygenic shoulder around a sparse causal); recommended for biobank-scale GWAS | Very small loci with truly sparse architecture (over-conservative); slower convergence |
| FINEMAP (Benner 2016 Bioinformatics 32:1493) | Shotgun stochastic search over causal configurations | .z + .ld + .master files | Exact Bayes factors at small k; widely cited | Slow at L > 5; binary install only (christianbenner.com); same LD-mismatch fragility as SuSiE |
| CAVIAR (Hormozdiari 2014 Genetics 198:497) | Exhaustive enumeration up to k causals | (z, R) | Exact posterior at small k | Combinatorial explosion beyond k=6; legacy method largely superseded by SuSiE |
| DAP-G (Wen 2016 AJHG 98:1114) | Deterministic posterior approximation with adaptive scan | SBAMS format; TORUS for enrichment priors | Fast at QTL scale (whole-transcriptome); pairs with TORUS hierarchical priors | SBAMS format is awkward; less ubiquitous tooling |
| PAINTOR (Kichaev 2014 PLoS Genet 10:e1004722) | EM with binary functional annotations | (z, R, A) per locus | Locus-level functional priors; multi-trait variant | Single-trait mode often matched by PolyFun + SuSiE; slower than SuSiE |
| PolyFun + SuSiE/FINEMAP (Weissbrod 2020 Nat Genet 52:1355) | Stratified LDSC genome-wide -> per-SNP prior_weights | GWAS sumstats + pre-baked baseline-LF | Most powerful single-trait functional prior; >20% more high-PIP (PIP>0.95) variants in simulations, >32% in real UK Biobank traits (Weissbrod 2020) | Requires matched-ancestry baseline-LF; runs in two stages |
| SuSiEx (Yuan 2024 Nat Genet 56:1841) | Joint cross-ancestry SuSiE; shared causal, population-specific LD | Per-pop sumstats + per-pop LD reference | Smaller credible sets than per-ancestry meta or marginal fine-mapping; principled when causal variants are shared | Trans-ethnic heterogeneity violated (population-specific causals); ancestry must be cleanly assigned |
| MultiSuSiE (Rossen 2025 Nat Genet) | Cross-ancestry SuSiE variant; flexible heterogeneity | Per-pop sumstats + per-pop LD | Similar to SuSiEx; alternative implementation | Same as SuSiEx; newer, less battle-tested |
| FOCUS / MA-FOCUS (Mancuso 2019 Nat Genet 51:675) | Probabilistic TWAS fine-mapping over gene models | TWAS Z-scores + gene LD (predicted expression) | Identifies likely causal gene among co-regulated TWAS hits; cross-ancestry MA-FOCUS variant | Requires pre-computed expression weights (e.g., FUSION/PrediXcan); gene-level rather than variant-level inference |

Methodology evolves; verify the latest susieR vignette and the SuSiE-inf paper before locking on a single method. Wang Lab maintains susieR; the IBSS algorithm is stable but argument semantics (e.g., `prior_weights` vs `prior_variance`) have changed across versions.

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| SuSiE finds 3 credible sets, FINEMAP finds 1 | FINEMAP's stochastic search did not converge OR SuSiE absorbed background into spurious sets | Increase FINEMAP `--n-iterations`; check SuSiE purity (sets with purity < 0.5 are spurious) |
| SuSiE PIPs much sharper than FINEMAP | susie_rss assumes single residual variance; FINEMAP marginalizes over noise | Both can be correct; report the intersection of high-PIP variants from both as primary candidates |
| PolyFun + SuSiE collapses 10-variant credible set to 1 | Functional priors are doing real work (coding variant in set) | Verify with `prior_weights` plot; if priors are coding-specific the result is interpretable |
| SuSiEx credible set excludes the EUR top-PIP variant | EUR signal is tag, true causal shared across ancestries lies elsewhere | Trust SuSiEx if both populations have well-powered GWAS; verify with conditional analysis |
| HLA gives 50-variant credible set in every method | HLA LD structure cannot be fine-mapped by linear methods | Use HLA-specific imputation (HIBAG, SNP2HLA) and haplotype-level analysis |
