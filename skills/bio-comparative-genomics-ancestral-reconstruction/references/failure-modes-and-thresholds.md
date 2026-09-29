## Per-Tool Failure Modes

### Long-branch attraction (LBA) at deep nodes

**Trigger:** Tree with two long terminal branches separated by a short internal branch; mixed amino-acid compositions across taxa.

**Mechanism:** Site-homogeneous models (LG, WAG, JTT) assume constant amino-acid equilibrium frequencies across the tree. When real compositions differ (e.g. thermophiles vs mesophiles), models underestimate the probability of convergent substitutions at compositionally-constrained sites, producing apparent shared derived states between long branches (Szánthó 2023 Syst Biol 72:767). The reconstructed ancestor at the deep node is biased toward whichever long-branch composition the model favors.

**Symptom:** Bootstrap support at the contested node remains high under site-homogeneous models but collapses under CAT-GTR / CAT-PMSF. Posterior predictive checks for compositional homogeneity reject the model (Foster 2004 Syst Biol 53:485).

**Fix:** Move to PhyloBayes-MPI with CAT-GTR or IQ-TREE2 with CAT-PMSF (`-m LG+C60+F+R` then `-ft <tree>` for posterior mean site frequencies). For ASR specifically, use ancestral reconstruction only when the model passes compositional adequacy. Slow-fast site removal (recoded amino acids; Susko & Roger 2007 MBE 24:2139) is an alternative.

### Epistasis breaking site-independent ASR

**Trigger:** Multiple sites in the same protein evolve under coupled constraints (compensatory pairs in RNA secondary structure; buried-residue covariance; allosteric networks).

**Mechanism:** ML/Bayesian ASR assumes sites are independent given the tree; the joint ancestral sequence is the product of per-site posteriors. Real proteins evolve through compensatory substitutions where a destabilizing mutation at site i is compensated by a mutation at site j (Pollock 2012 PNAS 109:E1352; Shah 2015 PNAS 112:E3226). The ML ancestral sequence can contain a never-tested combination of states.

**Symptom:** The reconstructed protein fails to fold or is non-functional when expressed; positions flagged ambiguous (P < 0.9) are non-random and cluster in 3D space when mapped to structure.

**Fix:** Use GRASP indel-aware reconstruction; design 4-8 alternative constructs varying ambiguous positions (P < 0.9), prioritizing residues that are structurally coupled to high-confidence ML states; experimentally test each construct; report the range of functional reconstructions, not the single ML sequence. Hochberg & Thornton 2017 Annu Rev Biophys 46:247 review epistasis strategies.

### Root placement error cascading to deep ancestors

**Trigger:** Trees rooted by midpoint, outgroup with very long branch, or `--prefix` auto-root.

**Mechanism:** Marginal ASR posteriors at internal nodes depend on the root's position because the root defines the time direction of substitution. A wrong root flips state inferences for deep nodes (especially when ancestral state is asymmetric, e.g. presence -> absence is more common than reverse).

**Symptom:** Re-rooting the tree changes the inferred ancestral state at the deepest node by > 0.2 posterior probability; STRIDE rooting (Emms 2017 MBE 34:3267) disagrees with outgroup rooting.

**Fix:** Run ASR over a set of candidate roots; report robust nodes (state invariant) and root-sensitive nodes separately. For phylogenomic-scale data, use STRIDE / MAD rooting (Tria 2017 Nat Eco Evo 1:0193) or ALE-rooting (Williams 2017 PNAS 114:E4602) and document the rooting strategy.

### BiSSE false-positive in state-dependent diversification

**Trigger:** Testing whether a discrete trait influences speciation/extinction using BiSSE (Maddison 2007 Syst Biol 56:701).

**Mechanism:** BiSSE attributes ALL rate variation to the focal trait. When real diversification heterogeneity is caused by a hidden character correlated with the focal trait, BiSSE reports a spurious significant association (Rabosky-Goldberg 2015 Syst Biol 64:340 -- ~40% Type-I rate at moderate trees).

**Symptom:** BiSSE LRT highly significant but biological mechanism unclear; HiSSE rejects BiSSE in favor of a hidden-state model with the focal trait neutral.

**Fix:** Run HiSSE as the required null model (Beaulieu & O'Meara 2016 Syst Biol 65:583). Report BiSSE only if HiSSE-null is rejected. For traits with deep clade structure, use FiSSE (Rabosky & Goldberg 2017 Evolution 71:1432) which is robust to model misspecification by design.

### Parsimony vs ML on asymmetric rates (Felsenstein bias)

**Trigger:** Trait has strongly asymmetric forward vs reverse rates (e.g. gene loss > gene gain).

**Mechanism:** Parsimony minimizes total changes, implicitly assuming symmetric rates. ML/Bayesian methods estimate the rate matrix from the data and reconstruct accordingly. Under strong asymmetry, parsimony over-reconstructs the rarer state at ancestors (Cunningham 1999 Syst Biol 48:665).

**Symptom:** Parsimony and ML reconstructions disagree at > 30% of nodes; ML rates fit AIC-better with ARD (all-rates-different) than ER (equal-rates).

**Fix:** Always run ER vs SYM vs ARD model comparison via `ape::ace(model='...')` AIC; use ARD when asymmetry is supported. For gene-content evolution, Dollo parsimony (gain rare, loss common) is often the better prior than equal-rates ML.

### Continuous-trait BM-only model with non-BM evolution

**Trigger:** Fitting `phytools::fastAnc()` (which assumes BM) to a trait with strong directional or stabilizing selection.

**Mechanism:** fastAnc returns the BM-MLE ancestral state, which is a weighted mean of descendant values with weights from the BM covariance matrix. If the trait evolved under OU (stabilizing), real ancestor values were closer to the optimum than fastAnc returns; if under EB (early burst), real ancestors were more variable than fastAnc returns.

**Symptom:** BM model fits with `geiger::fitContinuous(model='BM')` give AIC > 4 above OU or EB; phylogenetic signal Pagel's lambda < 0.5; Blomberg's K significantly < 1 (Blomberg 2003 Evolution 57:717).

**Fix:** Run `fitContinuous` with multiple models (BM/OU/EB/lambda/kappa/delta); use the best AIC model's ancestral reconstruction. For OU, use `OUwie::ace()`; for regime shifts, `bayou::bayou.mcmc`. Always report Pagel's lambda alongside ancestral estimates as a phylogenetic-signal indicator. Boettiger 2012 Evolution 66:2240 and Cooper 2016 Biol J Linn Soc 118:64 detail model-adequacy testing.

### Alignment error propagating to ASR

**Trigger:** Using `MUSCLE` or default MAFFT alignment on highly diverged sequences (< 30% identity).

**Mechanism:** Misaligned columns place non-homologous residues into the same column. ML ASR treats those residues as states of the same character, producing impossible ancestral inferences (Vialle 2018 MBE 35:1783).

**Symptom:** Ambiguous regions of the alignment correspond to low-confidence ASR sites; gappy columns dominate the low-confidence set; alignment scoring (TCS, Guidance2) marks the same regions as poorly aligned.

**Fix:** Filter alignment with HmmCleaner (Di Franco 2019 BMC Evol Biol 19:21) or PREQUAL (Whelan et al 2018 Bioinformatics 34:3929) before ASR. Segment-level filtering outperforms block-level filtering (Gblocks, trimAl) for downstream evolutionary inference. For ASR specifically, mask ambiguous columns (treat as missing) rather than removing them, to preserve coordinates.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| ASR site high confidence | posterior >= 0.95 | Standard convention (Yang 1995); above this treat state as fixed |
| ASR site moderate confidence | 0.80 <= posterior < 0.95 | Worth alternative-construct testing in resurrection studies |
| ASR site ambiguous | posterior < 0.80 | Design alternative constructs; cluster against structure |
| Pagel's lambda interpretation | lambda > 0.7 strong signal; 0.3-0.7 moderate; < 0.3 weak (ad-hoc operational convention; Pagel 1999 introduced lambda but did not prescribe these cutoffs) | Pagel 1999 Nature 401:877 (method); community convention (thresholds) |
| Blomberg K interpretation | K > 1 conserved; K = 1 BM; K < 1 weak signal | Blomberg 2003 Evolution 57:717 |
| Bootstrap support for ancestral clade | >= 70% before trusting the state at that node | Standard; below this, root-sensitivity tests required |
| MCMC ESS for Bayesian ASR | ESS >= 200 per parameter; ASRV at least 200 | RevBayes / Tracer convention; Lakner 2008 Syst Biol 57:86 |
| Stochastic mapping nsim | >= 1000 simulations per tree | Operational convention (Bollback 2006 BMC Bioinf 7:88 SIMMAP method); for asymmetric rates raise to >= 5000 |
| Tree depth limit for protein ASR | dS / branch length < 1.0 at deepest node | Above this, signal saturated; Yang 2007 PAML manual |
| Codon ASR minimum sequences | >= 8 with sufficient divergence (~0.5 substitutions/site total) | Operational convention; below this, codon-model parameters poorly constrained |
| Minimum taxa for binary trait ER vs ARD AIC | >= 20 tips; below this, rates often unidentifiable | Beaulieu 2016 |
| GRASP indel posterior threshold | >= 0.8 to call indel present at node | GRASP documentation; below this, both states tested experimentally |
| OUwie regime requires | >= 10 tips per regime | Beaulieu 2012 Evolution 66:2369; below this, optima unidentifiable |
| HMM rate categories | start with rate.cat=2; AIC compare against 1 | Boyko & Beaulieu 2021 MEE 12:468 |

