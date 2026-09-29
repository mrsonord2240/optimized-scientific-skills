## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Parsimony and ML disagree at > 30% nodes | Rate asymmetry, model misspecification, or weak signal may contribute | Compare supported transition models and check long-branch/root sensitivity; report the competing node states and uncertainty rather than selecting by a fixed disagreement threshold |
| BiSSE highly significant, HiSSE neutral | Hidden trait drives diversification | Report HiSSE as primary (Rabosky-Goldberg 2015) |
| Codeml marginal and joint disagree at deep node | Strong epistasis or model misspecification | Test alternative constructs; consider GRASP and BAli-Phy joint inference |
| Site-homogeneous LG agrees with PhyloBayes CAT-GTR at all nodes | Agreement across these fits on the tested data | Report both fits and their assumptions; agreement does not establish absence of compositional heterogeneity or general model adequacy |
| Site-homogeneous LG and CAT-GTR disagree at deep node | Composition, site heterogeneity, rooting, or other model differences may contribute | Retain both estimates, inspect model fit and root sensitivity, and report the node as model-sensitive; CAT-PMSF may be a follow-up comparison where appropriate |
| Marginal and stochastic-map states disagree | Rate matrix, root prior, conditioning, or finite mapping variation may contribute | Report both summaries, document the root prior (including `pi='fitzjohn'` where used), and assess mapping stability; neither method is automatically preferred |
| BM ASR vs OU ASR disagree on direction of trait change at root | The ancestral estimate is sensitive to the model assumptions | Compare model adequacy, data support, regime specification, and uncertainty. AIC support alone does not establish that OU is correct; report both estimates and the sensitivity |
| Reconstruction at deepest node flips under re-rooting | Insufficient outgroup support; LBA | Run STRIDE / MAD rooting; report deep-node state with uncertainty |
| GRASP indel and PAML "missing-data" reconstructions disagree | The methods encode gaps and indels differently | Describe each method's gap treatment and retain the difference as uncertainty; choose based on the scientific question and validate alternative constructs rather than treating either output as ground truth |

**Operational rule for publication:** Compare scientifically plausible models (for example, ER/SYM/ARD for discrete traits, BM/OU/EB for continuous traits, and site-homogeneous with a site-heterogeneous comparison for deep sequence ASR) when data and validated tools support them. Report the assumptions and uncertainty for each result. If reasonable models differ, retain competing estimates and explicitly flag model-sensitive nodes; reconstructed states are estimates, not observed data.

## Cohort Gotchas

- **Polyploid species in continuous-trait analyses:** body size, genome size, gene count are confounded with ploidy; assign subgenomes (see [[whole-genome-duplication]]) and treat as separate tips, or use multilabel-tree methods.
- **Hybrid taxa break the bifurcating-tree assumption:** ape and phytools assume strict bifurcations; for hybrids, use phylogenetic networks (phangorn, RevBayes admixture) or remove hybrid tips before ASR.
- **Tip-dated trees from molecular clock require careful root prior:** RevBayes / BEAST2 with calibrated tip dates produce trees in absolute time; PAML/IQ-TREE work in relative substitutions. Match the time unit when integrating across tools.
- **OrthoFinder species trees from gene-tree summary (STAG):** branch lengths are coalescent-units when used for ASR; convert to substitutions via concatenated alignment if downstream tools require it.

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why marginal not joint reconstruction?" | Marginal exposes per-site uncertainty necessary for resurrection construct design; joint is internally consistent but hides ambiguity (Pupko 2000 MBE 17:890) |
| "How was epistasis handled?" | Designed and tested N alternative constructs at ambiguous (P < 0.8) sites; report functional range, not just ML sequence (Hochberg & Thornton 2017) |
| "Why these models?" | AIC compared ER/SYM/ARD for discrete; BM/OU/EB/lambda for continuous; site-homogeneous + CAT-PMSF for deep sequence ASR; reported ancestral state only at model-invariant nodes |
| "Phylogenetic signal?" | Pagel's lambda = X; Blomberg's K = Y; signal supports tree-based ASR (or: signal weak, ASR exploratory only) |
| "Effect of rooting?" | Reconstructed under multiple rootings; state at root invariant across STRIDE / MAD / outgroup, or explicit caveat for root-sensitive nodes |
| "Multiple-testing across nodes?" | Ancestral states are estimated quantities, not tested hypotheses; report posteriors per node, no multiplicity correction needed |
| "Why not BiSSE for trait-diversification correlation?" | HiSSE replaces BiSSE per Rabosky-Goldberg 2015; BiSSE Type-I rate ~40% on simulated data |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `codeml` rst file missing posteriors section | `RateAncestor` set to 0 or not set | Set `RateAncestor=1` in control file |
| IQ-TREE `--ancestral` empty `.state` file | Tree not rooted | Specify outgroup with `-o`; IQ-TREE requires rooting for ancestral output |
| `make.simmap` chains never mix | Asymmetric rate matrix with sparse data | Increase nsim to 5000+; switch root prior to `pi='fitzjohn'` |
| `ape::ace` fails with `NA/NaN/Inf in foreign function call` | Polytomies in tree | `multi2di()` to resolve to bifurcating; or use `ape::ace(method='ML')` with phytools |
| GRASP runs but produces empty asr.fasta | Alignment includes stop codons or X characters | Strip non-canonical residues; check that the alignment passes Bio.SeqIO validation |
| `fitContinuous(model='OU')` returns lambda = 0 | OU collapsed to white noise (no signal) | Trait variance unexplained by tree; reconsider whether continuous-trait ASR is meaningful |
| Stochastic mapping returns all-or-nothing at deep nodes | Insufficient signal; pi='equal' default | Use `pi='fitzjohn'`; check that taxa span both states |
| HiSSE convergence failures | `bound.par` too restrictive; hessian singular | Use `starting.vals=NULL`, restart with `output.liks=TRUE`; switch to `BiSSE-ness` (Magnuson-Ford & Otto 2012 Am Nat 180:225) as fallback |
| codeml omega = 0 or 999 at branch | Saturation or alignment artifact | Increase model complexity (M0 -> M3); check dS at branch; if dS > 3, ASR unreliable on that branch |
| Ancestral genome content reconstruction inflated | Assembly fragmentation produces false absences | Use BUSCO-completeness-corrected presence/absence; or run [[gene-tree-species-tree-reconciliation]] which handles loss-vs-missing |

## Tool Installation Notes

```bash
# CLI
conda install -c bioconda paml iqtree
# RevBayes via source build (https://revbayes.github.io/download) or homebrew on macOS
# GRASP from https://github.com/bodenlab/GRASP (Java)
# FastML web at fastml.tau.ac.il or CLI binary

# R
install.packages(c('ape', 'phytools', 'geiger', 'phangorn', 'OUwie', 'bayou', 'RERconverge'))
# Tested source for the shipped corHMM example: CRAN 2.8
remotes::install_version('corHMM', version='2.8', repos='https://cran.r-project.org')
remotes::install_github('thej022214/hisse')

# Python
pip install biopython ete4
```

For protein resurrection, consider GRASP when indel-aware reconstruction is important. For selection-context codon ASR, PAML codeml is an established option. For trait macroevolution, choose among phytools, corHMM, and OUwie according to the trait type, model assumptions, data support, and sensitivity results; hidden-rate corHMM models are useful when rate heterogeneity is plausible and estimable. Use Bayesian methods such as RevBayes when the question requires a full posterior or explicit model uncertainty and the analysis can support the added computation.

## References

- Yang Z et al 1995 Genetics 141:1641 (marginal ASR likelihood framework)
- Pupko T et al 2000 MBE 17:890 (joint ASR efficient algorithm)
- Nielsen R 2002 Syst Biol 51:729 (stochastic mapping)
- Huelsenbeck JP et al 2003 Syst Biol 52:131 (Bayesian stochastic mapping)
- Felsenstein J 1985 Am Nat 125:1 (phylogenetic independent contrasts)
- Felsenstein J 2012 Am Nat 179:145 (threshold model)
- Pagel M 1999 Nature 401:877 (lambda phylogenetic signal)
- Blomberg SP et al 2003 Evolution 57:717 (K statistic)
- Beaulieu JM et al 2013 Syst Biol 62:725 (hidden Markov state-rate decoupling)
- Beaulieu JM & O'Meara BC 2016 Syst Biol 65:583 (HiSSE)
- Rabosky DL & Goldberg EE 2015 Syst Biol 64:340 (BiSSE Type-I rates)
- Boyko JD & Beaulieu JM 2021 MEE 12:468 (generalized HMM corHMM)
- Bollback JP 2006 BMC Bioinf 7:88 (SIMMAP)
- Pollock DD et al 2012 PNAS 109:E1352 (compensatory epistasis)
- Shah P et al 2015 PNAS 112:E3226 (contingency and entrenchment epistasis)
- Hochberg GKA & Thornton JW 2017 Annu Rev Biophys 46:247 (ASR for protein resurrection)
- Foley G et al 2022 PLoS Comp Biol 18:e1010633 (GRASP indel-aware ASR)
- Szánthó LL et al 2023 Syst Biol 72(4):767-780 (compositional LBA; CAT-PMSF) -- DOI 10.1093/sysbio/syad013
- Boettiger C et al 2012 Evolution 66:2240 (model adequacy for continuous-trait macroevolution)
- Cooper N et al 2016 Biol J Linn Soc 118:64 (cautionary note OU)
- Cunningham CW 1999 Syst Biol 48:665 (asymmetric rate parsimony bias)
- Felsenstein J 1978 Syst Zool 27:401 (long branch attraction)
- Maddison WP et al 2007 Syst Biol 56:701 (BiSSE)
- Beaulieu JM et al 2012 Evolution 66:2369 (OUwie)
- Uyeda JC & Harmon LJ 2014 Syst Biol 63:902 (bayou)
- FitzJohn RG 2009 Syst Biol 58:595 (root prior)
- Whelan S et al 2018 Bioinformatics 34:3929 (PREQUAL)
- Di Franco A et al 2019 BMC Evol Biol 19:21 (HmmCleaner)
- Emms DM & Kelly S 2017 MBE 34:3267 (STRIDE rooting)
- Tria FDK et al 2017 Nat Eco Evo 1:0193 (MAD rooting)
- Williams TA et al 2017 PNAS 114:E4602 (ALE-rooting of deep phylogenies)
- Hu Z et al 2019 MBE 36:1086 (PhyloAcc)
- Fukushima K & Pollock DD 2023 Nat Eco Evo 7:155 (CSUBST)
- Muffato M et al 2023 Nat Eco Evo 7:355 (AGORA)

## Related Skills

- comparative-genomics/positive-selection - Branch- and site-level selection inference on ancestral branches
- comparative-genomics/ortholog-inference - Define orthogroups whose alignments feed ASR
- comparative-genomics/gene-tree-species-tree-reconciliation - DTL-aware ancestral gene-content inference; root inference via ALE
- comparative-genomics/whole-genome-duplication - Ks-dating provides time scale for ancestral state inference
- comparative-genomics/comparative-annotation-projection - Project ancestral CDS to descendants for validation
- phylogenetics/modern-tree-inference - Generate rooted ML/Bayesian trees as ASR scaffold
- phylogenetics/bayesian-inference - RevBayes / MrBayes priors for Bayesian ASR
- phylogenetics/divergence-dating - Time-calibrated trees as input for absolute-time ASR
- alignment/multiple-alignment - PRANK / MACSE indel-aware alignment before sequence ASR
- alignment/alignment-trimming - PREQUAL / HmmCleaner filtering before ASR
