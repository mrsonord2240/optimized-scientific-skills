## PAML codeml Ancestral Reconstruction

**Goal:** Reconstruct ancestral codon or protein states at internal nodes using ML under a stationary codon/protein model.

**Approach:** Build a codon-aware alignment (PRANK or MACSE) -> infer rooted ML tree with the same model intended for ASR -> create codeml control file with `RateAncestor=1` -> run codeml -> parse the current matrix-oriented `rst` marginal section for each node's best state and its probability. This value is not a complete distribution across alternative states; retain the limitation when reporting uncertainty. The parser checks the node and site records and raises an error if they are absent or malformed.

[Run or adapt `codeml_asr.py`](../scripts/codeml_asr.py) for the reusable control-file and posterior-parsing example.


## IQ-TREE2 Marginal Reconstruction

**Goal:** Faster ASR with native model selection, for protein/DNA alignments where PAML is too slow.

**Approach:** Run IQ-TREE with `-m TEST` to pick model -> `--ancestral` produces `.state` table of per-site posteriors; output rooted by `-o <outgroup>`.

[Run or adapt `iqtree_ancestral.sh`](../scripts/iqtree_ancestral.sh) for the IQ-TREE2 command example.


[Run or adapt `iqtree_state.py`](../scripts/iqtree_state.py) for the `.state` parser example.


## Stochastic Mapping (phytools::make.simmap)

**Goal:** Sample full character histories along branches for a discrete trait; quantify ancestral state probabilities with proper uncertainty.

**Approach:** Fit Mk rate matrix -> sample `nsim` simulated maps under the posterior -> summarize state probabilities per node and transitions per branch.

[Run or adapt `stochastic_mapping.R`](../scripts/stochastic_mapping.R) for the Mk-model comparison, stochastic-mapping, and hidden-rate example.

Run it with an explicit seed, for example `Rscript stochastic_mapping.R --seed=20260928`. The script rejects mismatched tree/data taxa and missing, blank, or single-state trait inputs before fitting, and prints seed metadata with the result.


`pi='fitzjohn'` uses the Fitzjohn 2009 Syst Biol 58:595 root prior (root-state prior weighted by the likelihood of the observed tip data given each root state), preferable to `pi='estimated'` which fixes the prior to the estimated equilibrium distribution and can over-fit, or `pi='equal'` which can bias toward the rarer state when data are asymmetric.

## Continuous-Trait ASR with Model Adequacy

**Goal:** Reconstruct ancestral values for a continuous trait while honestly reporting which model class the data support.

**Approach:** Fit and compare BM/OU/EB/lambda/kappa/delta, then reconstruct only under a method that applies the selected fitted model. The packaged example supports BM `fastAnc` confidence intervals. OUwie 3.0.3 accepts a fitted OUwie object and `knowledge=TRUE`, but its ancestral estimates have no uncertainty and are for visualization/model intuition; the script labels them exploratory. For EB/lambda/kappa/delta winners, the example stops and requests a method that reconstructs under the fitted transformation instead of substituting original-tree BM output. Report Pagel's lambda and Blomberg's K as signal summaries, not proof of model adequacy.

[Run or adapt `continuous_trait_asr.R`](../scripts/continuous_trait_asr.R) for the BM/OU/EB model-comparison and reconstruction example.


Report Pagel's lambda alongside ancestral values as a descriptive signal statistic. Any numeric interpretation threshold is an operational convention; it does not by itself establish BM-like evolution, identify ecological causes, or show that an ASR is well constrained.

## GRASP Indel-Aware Sequence ASR (Protein Resurrection Workflow)

**Goal:** Reconstruct ancestral protein sequence INCLUDING indel states, for experimental resurrection.

**Approach:** GRASP (Foley 2022 PLoS Comp Biol 18:e1010633) uses a partial-order alignment graph to model indels probabilistically; outputs alternative reconstructions and explicit indel uncertainty.

[Run or adapt `grasp_asr.sh`](../scripts/grasp_asr.sh) for the GRASP command and output-surface example.


**Construct-design protocol:** (1) Take ML sequence as primary construct. (2) For each site with max posterior < 0.8 and a second state with posterior > 0.2, build a single-mutant alternative. (3) For each indel block with posterior < 0.8, build a present and an absent alternative. (4) Order constructs by structural compactness (avoid surface loops first). Typical batch: 4-8 constructs per ancestral node. Hochberg & Thornton 2017 Annu Rev Biophys 46:247 review the operational pipeline.

