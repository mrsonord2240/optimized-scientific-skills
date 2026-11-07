# PP.H4 Threshold Framework
| Threshold | Use case | Source |
|-----------|----------|--------|
| 0.5 - 0.7 | Suggestive / pilot / hypothesis-generating | Giambartolomei 2014 original |
| >= 0.7 | Triangulation tier for TWAS / cis-MR / effector-gene cross-evidence | Open Targets Genetics common practice; cross-reference downstream skills |
| >= 0.75 | Open Targets Platform / eQTL Catalogue / FinnGen default screening threshold | Open Targets Genetics docs; Mountjoy 2021 Nat Genet 53:1527 |
| >= 0.80 | Most published colocalizations / standard publication tier | Wallace 2020 PLoS Genet 16:e1008720 |
| >= 0.90 | Stringent clinical / therapeutic-target prioritization | Reserved for high-confidence claims |
| >= 0.95 | Industry / regulatory drug-target submission grade | Internal pharma default |
| PP.H3 >= 0.80 | Confident distinct causal variants (negative coloc result) | Standard |
| PP.H4 / (PP.H3 + PP.H4) >= 0.9 | Conditional probability framing (some pipelines) | Foley 2021 |

**Operational rule:** Three operational tiers map onto the most common downstream uses: (a) **>= 0.7** when PP.H4 is one of several lines of triangulating evidence (TWAS + coloc, cis-MR + coloc, effector-gene multi-evidence) -- this is the threshold downstream skills (causal-genomics/transcriptome-wide-association, causal-genomics/mendelian-randomization cis-MR, causal-genomics/effector-gene-prioritization, causal-genomics/proteome-mr-drug-target) require; (b) **>= 0.8** for standard peer-reviewed publication as a stand-alone coloc claim (Wallace 2020); (c) **>= 0.95** for industry / clinical drug-target submission. Open Targets and FinnGen pipelines screen at >= 0.75 but downstream publication-grade coloc claims should clear >= 0.8 and triangulation claims >= 0.7. ALWAYS report PP.H3 alongside PP.H4 -- a locus with PP.H4 = 0.6, PP.H3 = 0.3 is qualitatively different from PP.H4 = 0.6, PP.H3 = 0.05 (the former is real ambiguity over single vs distinct causal; the latter is underpowered evidence). Run `coloc::sensitivity()` and report the p12 range over which PP.H4 stays above the threshold.

## Default Priors and the p12 Sensitivity Question

| Prior | Default | Interpretation | When to change |
|-------|---------|----------------|----------------|
| p1 | 1e-4 | Prob a random SNP is associated with trait 1 | Rarely changed |
| p2 | 1e-4 | Prob a random SNP is associated with trait 2 | Rarely changed |
| p12 | 1e-5 | Prob a random SNP is associated with both traits | Lower (5e-6 or 1e-6) for trans-eQTL or unrelated trait pairs; raise (5e-5) only with strong prior, e.g. molecular QTL in the same tissue as causal cell type |

The p12/p1 ratio (= 0.1 under defaults) is the prior odds of colocalization given a trait-1 association. Wallace 2020 (PLoS Genet 16:e1008720) showed default p12 = 1e-5 is too liberal for many real-world settings and recommended sensitivity analysis as standard practice. Pullin & Wallace 2025 (PLoS Genet 21:e1011697) extended this with variant-specific priors weighted by functional annotation.

### p12 Sensitivity Grid

| p12 grid point | Use case | Reporting rule |
|----------------|----------|-----------------|
| 1e-4 | Suggestive only / EUR cis-eQTL relaxed | PP.H4 here cannot support a publication claim |
| 1e-5 | Default for most cis-eQTL <-> GWAS pairs | Standard |
| 5e-6 | Conservative cis; default for trans-eQTL coloc | Recommended publication baseline |
| 1e-6 | Very conservative; trans coloc with weak prior | Required for cross-trait genome-rare coloc |

**Operational rule:** Require PP.H4 to remain above threshold across at least 3 adjacent grid points; report the lowest p12 at which PP.H4 >= 0.75. Use `coloc::sensitivity(result, rule = 'H4 > 0.75')` for the diagnostic plot.

Required reporting: PP.H4 at default priors + p12 range over which PP.H4 stays above threshold.
