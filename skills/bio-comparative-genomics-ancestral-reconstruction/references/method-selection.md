## Version Compatibility

Version/source support for executable examples: PAML 4.10.10, IQ-TREE 2.4.0, official GRASP CLI dated 21-Mar-2024, R 4.4.3, ape 5.8-1, phytools 2.5-2, geiger 2.0.12, OUwie 3.0.3, and corHMM 2.8 installed from the CRAN source archive. The corHMM stochastic-mapping example is validated on CRAN 2.8; do not infer compatibility with other releases from this result.

The corHMM GitHub repository revision inspected on 2026-09-28 is commit
`3ae10b22cfb519245b75777342812a0023ac0d43` (DESCRIPTION version 2.10.5).
That revision adds RTMB and was not installed or executed in the prepared
environment; it is recorded for provenance, not presented as a supported or
tested interface. Install the tested CRAN source with
`remotes::install_version("corHMM", version = "2.8", repos = "https://cran.r-project.org")`.
Validate the package version and the exact example call again if changing the
source or version.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('corHMM')` then `?ancRECON`, `?make.simmap`, `?ace`
- CLI: `codeml` (no `--version`; check `paml -h` or examine `Phylip.tre` example), `iqtree2 --version`, `revbayes --version`
- Python: `pip show biopython`; check `Bio.Phylo.PAML.codeml` API

If code throws `AttributeError`, `ImportError`, missing slot errors on R S4 objects, or PAML `mlc` parsing failures, introspect the installed package (`?` in R, `help()` in Python) and adapt the example rather than retrying. PAML output formats are stable across 4.9 -> 4.10; IQ-TREE's `--ancestral` flag replaced `-asr` in v2.0+.

# Ancestral State Reconstruction

**"What did this gene / trait / genome look like at an internal node?"** -> Choose the reconstruction framework that matches the data class (sequence / discrete trait / continuous trait / gene content) and the inference question (point estimate vs full posterior; marginal vs joint vs scaled-conditional). The single most common mistake is reconstructing under a site-independent or trait-stationary model when the underlying biology demands a hidden-rate or epistatic model -- the resulting "ancestor" is mathematically optimal under the wrong model and is silently wrong (Beaulieu & O'Meara 2016 Syst Biol 65:583; Boyko & Beaulieu 2021 MEE 12:468).

- Sequence ASR (protein resurrection): PAML codeml `RateAncestor=1`; IQ-TREE2 `--ancestral`; GRASP (graph-based, handles indels); FastML (Bayesian)
- Discrete traits: R `ape::ace(type='discrete')`; `corHMM::corHMM()` (rate categories); `phytools::make.simmap()` stochastic mapping; BayesTraits MultiState/Discrete
- Continuous traits: `phytools::fastAnc()`; `phytools::contMap()`; `geiger::fitContinuous(model='BM'|'OU'|'EB')`; RPANDA `fit_t_env()`
- Ancestral gene content (presence/absence): `ape::ace(type='discrete', model='ARD')`; Dollo parsimony in phangorn; ALE/GeneRax for full DTL (see [[gene-tree-species-tree-reconciliation]])

## Algorithmic Taxonomy

| Framework | Data class | Inference | Strength | Fails when |
|-----------|------------|-----------|----------|------------|
| ML marginal (codeml RateAncestor=1; IQ-TREE --ancestral) | Sequence / discrete | Site-by-site MAP + posterior | Per-site uncertainty; fast (Yang 1995 Genetics 141:1641; Pupko 2000 MBE 17:890) | Strong epistasis; site-independent assumption violated |
| ML joint (Pupko 2000 MBE 17:890; IQ-TREE marginal+joint output) | Sequence / discrete | Single most-likely joint history across all nodes | Internally consistent ancestral sequence | Loses per-site uncertainty; epistasis hidden |
| Stochastic mapping (Nielsen 2002 Syst Biol 51:729; Huelsenbeck 2003 Syst Biol 52:131; phytools::make.simmap) | Discrete | Full posterior over character histories along branches | Quantifies transition timing and rates per branch; supports posterior arithmetic | Long trees (mixing slow); rare-state biases |
| Bayesian MCMC (RevBayes, BayesTraits, MrBayes) | Sequence / discrete / continuous | Full posterior; supports model averaging | Honest uncertainty; rate-variable; hierarchical | Slow; convergence diagnostics required (ESS > 200) |
| Parsimony (Fitch 1971 Syst Zool 20:406; Dollo) | Discrete | MP states at nodes | Fast; assumption-light | Felsenstein-zone LBA artifact; biased toward fast change (Felsenstein 1978 Syst Zool 27:401) |
| Hidden Markov / hidden rates (corHMM; HiSSE; HMM) | Discrete | State + rate-class jointly | Captures rate heterogeneity across the tree; non-stationarity (Beaulieu 2013 Syst Biol 62:725) | Requires enough state changes to identify hidden rates |
| Threshold model (Felsenstein 2012 Am Nat 179:145; phytools::threshBayes) | Binary on continuous liability | MCMC on latent liabilities | Models polygenic / underlying-quantitative discrete traits | Slow MCMC; complex liability covariance |
| BM / OU / EB on continuous (geiger::fitContinuous) | Continuous | Phylogenetic regression on BM, OU mean-reverting, EB time-decay | Standard for body-size / niche-shape continuous traits | Model adequacy ignored (Boettiger 2012 Evolution 66:2240; Cooper 2016 Biol J Linn Soc 118:64) |
| Multi-rate BM / OUwie (Beaulieu 2012 Evolution 66:2369) | Continuous | Rate / optimum varies by clade or discrete regime | Models regime shifts; integrates with discrete trait history | Regime mismapping cascades to spurious rate differences |
| Phylogenetic generalized least squares (PGLS) | Continuous (multivariate) | Mean expected under BM; covariance from tree | Tests for correlation while controlling shared ancestry (Felsenstein 1985 Am Nat 125:1) | Strong evolutionary rate heterogeneity; non-BM trait |
| DTL reconciliation for gene content (ALE, GeneRax) | Gene tree / orthogroup | Ancestral gene presence + duplications/transfers/losses | Joint sequence + gene-content posterior | See [[gene-tree-species-tree-reconciliation]] |
| Indel-aware ASR (GRASP, FastML) | Sequence | Treats gaps as a separate process | Handles indel evolution explicitly; supports protein engineering | Slower; limited model families |

Methodology evolves; verify the latest `corHMM` / `phytools` vignettes before locking on a single approach. For continuous-trait macroevolution, consult Cooper 2016 model-adequacy reviews.

## Decision Tree by Experimental Scenario

| Scenario | Recommended method | Why |
|----------|---------------------|-----|
| Protein resurrection (~50-500 Myr divergences) | IQ-TREE2 `--ancestral` + GRASP indel reconstruction | Per-site marginal probabilities for alt-construct design; GRASP fixes indel ambiguity that PAML treats as missing data |
| Codon-level sequence ASR with selection inference | PAML codeml `RateAncestor=1`, model M0 (single omega), `seqtype=1` | Codon model native; integrates with branch reconstruction; produces `rst` with BEB-style site probs |
| Deep eukaryote / archaeal ASR (> 1 Bya) | Bayesian (RevBayes / PhyloBayes-MPI CAT-GTR) | Site-heterogeneous CAT model corrects compositional LBA (Szánthó 2023 Syst Biol 72:767); ML site-homogeneous models fail at this depth |
| Binary discrete trait with 5-30 taxa | `ape::ace(type='discrete', model='ARD')` + bootstrap | Standard for simple binary; ER/SYM/ARD model comparison via AIC |
| Binary discrete trait with 30+ taxa, suspected rate variation | Compare `corHMM(rate.cat=2)` with suitable Mk alternatives | Consider hidden rates when the biology and observed changes support estimable heterogeneity; a sensitivity result motivates model-adequacy checks but does not make a hidden-rate model mandatory |
| State-dependent diversification (correlation with speciation/extinction) | HiSSE (Beaulieu & O'Meara 2016) NOT BiSSE | BiSSE has catastrophic Type-I rate when rate heterogeneity is misattributed (Rabosky & Goldberg 2015 Syst Biol 64:340); HiSSE is the required null |
| Multi-state with phylogenetic uncertainty | `phytools::make.simmap(nsim=1000)` over a tree distribution | Marginalize over tree + state uncertainty; report posterior probabilities |
| Continuous trait, single regime | `phytools::fastAnc()` + `contMap` | Fast BM ML reconstruction; visual continuous reconstruction along branches |
| Continuous trait, suspected regime shifts | OUwie or `bayou` (Uyeda 2014 Syst Biol 63:902) | Multi-optimum OU models infer optimum shifts and their tree positions |
| Binary trait expected to be polygenic underlying | Threshold model `phytools::threshBayes` | Models latent liability properly; binary -> continuous bridge (Felsenstein 2012) |
| Ancestral gene family content | DTL reconciliation (ALE / GeneRax) | See [[gene-tree-species-tree-reconciliation]]; full posterior over D/T/L events |
| Ancestral genome architecture (gene order) | AGORA (Muffato 2023 Nat Eco Evo 7:355); DeCoSTAR | Joint reconciliation + adjacency posterior |
| Convergent rate shifts in noncoding | PhyloAcc (Hu 2019 MBE 36:1086); Thomas 2024 update | Bayesian Markov model on conserved noncoding elements |
| Convergent amino-acid substitutions | CSUBST (Fukushima & Pollock 2023 Nat Eco Evo 7:155) | Combinatorial substitution ratio omega_C; null-corrected |
| Categorical trait correlated rates | RERconverge (Kowalczyk 2019; Redlich 2024 MBE 41:msae210) | Relative evolutionary rates linked to a binary or categorical phenotype |

