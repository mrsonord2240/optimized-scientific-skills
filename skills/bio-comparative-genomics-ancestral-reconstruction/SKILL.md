---
name: bio-comparative-genomics-ancestral-reconstruction
description: Reconstruct ancestral states at internal phylogenetic nodes for sequences (PAML codeml, IQ-TREE --ancestral, GRASP, FastML), discrete traits (corHMM hidden-rate Markov, ape::ace, phytools::make.simmap stochastic mapping, BayesTraits), and continuous traits (phytools::fastAnc, geiger Brownian/OU, RPANDA). Use when designing constructs for ancestral protein resurrection, tracing trait evolution along a tree, performing stochastic character mapping, testing models of trait evolution (BM vs OU vs EB), inferring ancestral genome content via Dollo or DTL reconciliation, or quantifying ancestral-state uncertainty for downstream comparative analyses.
tool_type: mixed
primary_tool: PAML
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples target PAML 4.10.7+, IQ-TREE 2.3.6+, GRASP 2024+,
FastML 3.11+, RevBayes 1.2.4+, BayesTraits V4.1+, R 4.4+, ape 5.8+,
phytools 2.3+, corHMM 2.8 (CRAN; tested), geiger 2.0.11+, phangorn 2.12+,
RERconverge 0.3.0+, and Biopython 1.84+.

Before using a pattern, verify the installed version and inspect the relevant
CLI help or package documentation. If an API or output format differs, adapt
the example to the installed version rather than retrying unchanged. Read the
full provider version matrix in
[method-selection.md](references/method-selection.md) and symptom guidance in
[failure-modes-and-thresholds.md](references/failure-modes-and-thresholds.md).
For corHMM, this package's executable example is verified with CRAN 2.8 only;
see the source and support matrix in `method-selection.md`. A different source
or version requires independent interface validation before use.

# Ancestral State Reconstruction

Match the reconstruction framework to the data class and inference question.
Decide explicitly whether the goal is a point estimate or full posterior and
whether marginal, joint, or scaled-conditional states are required. Do not
interpret an ancestor reconstructed under a site-independent, stationary, or
Brownian model as biologically robust until that model and the root placement
have been checked.

## Core Workflow

1. Classify the input as sequence, discrete trait, continuous trait, or gene
   content. Name the internal node, state representation, and downstream use.
2. Validate the input before reconstruction: check taxon matching, missingness,
   alignment quality, tree topology and branch lengths, rooting, polytomies,
   trait coding, and enough variation to identify the selected model.
3. Select at least one primary method and a scientifically distinct comparison
   model. Use the scenario table in
   [method-selection.md](references/method-selection.md); do not choose a method
   solely because its software is available.
4. Record software versions, model, priors or rate matrix, root treatment,
   random seed, and all thresholds. Run the packaged example only after
   adapting paths and declared inputs.
5. Quantify uncertainty at the same resolution as the claim: per-site state
   posteriors, node-state probabilities, stochastic histories, confidence
   intervals, or posterior distributions. Preserve ambiguous states rather
   than silently collapsing them.
6. Test sensitivity to plausible models, alignments, roots, tree uncertainty,
   rate heterogeneity, and missing data. For deep sequence reconstruction,
   check compositional adequacy and site heterogeneity.
7. Reconcile disagreements using
   [reporting-and-troubleshooting.md](references/reporting-and-troubleshooting.md).
   Report model-sensitive nodes explicitly and downgrade single-model claims.
8. Deliver the reconstructed states together with inputs, exact commands or
   scripts, model-fit evidence, uncertainty, sensitivity results, and the
   limits on any resurrection or evolutionary interpretation.

## Choose the Surface

| Question | Primary surface | Required comparison or check |
|---|---|---|
| Protein or codon sequence | PAML codeml or IQ-TREE2; GRASP for indels | Alternative models, roots, and ambiguous-site constructs |
| Discrete trait | `ape::ace`, `phytools::make.simmap`, or corHMM | ER/SYM/ARD comparison and hidden-rate sensitivity |
| Continuous trait | `phytools::fastAnc` plus `geiger::fitContinuous` | BM/OU/EB comparison and phylogenetic signal |
| Full posterior | RevBayes or BayesTraits | Convergence and effective sample size |
| Gene content | Dollo, presence/absence ML, or DTL reconciliation | Missingness and assembly-completeness sensitivity |

Read [method-selection.md](references/method-selection.md) before selecting
between marginal and joint reconstruction, stochastic mapping, hidden-rate
models, threshold models, BM/OU/EB, or DTL reconciliation.

## Non-Negotiable Guardrails

- Do not trust a deep sequence ancestor when the alignment, root, or
  compositional model is unstable.
- Do not treat site-independent maximum-likelihood residues as an experimentally
  verified protein. Preserve posterior ambiguity and design alternative
  constructs for uncertain or structurally coupled positions.
- Do not use parsimony as the sole result when forward and reverse transition
  rates are asymmetric.
- Do not report a Brownian-motion continuous-trait reconstruction without
  comparing plausible alternatives and describing phylogenetic signal.
- Do not infer trait-dependent diversification from BiSSE alone; use a
  hidden-state null such as HiSSE as described in the routed guidance.
- Do not turn operational thresholds into universal biological constants.
  Report their source, purpose, and sensitivity.
- Do not report only the preferred model when reasonable models change the
  ancestral state. Flag those nodes as model-sensitive.

## Routed Resources

- Read [method-selection.md](references/method-selection.md) for the full
  algorithm taxonomy and scenario-specific decision table.
- Read
  [failure-modes-and-thresholds.md](references/failure-modes-and-thresholds.md)
  for long-branch attraction, epistasis, root errors, rate asymmetry,
  non-Brownian evolution, alignment errors, and quantitative thresholds.
- Read [method-workflows.md](references/method-workflows.md) for the PAML,
  IQ-TREE2, stochastic-mapping, continuous-trait, and GRASP workflows. Its code
  is packaged under `scripts/`.
- Read
  [reporting-and-troubleshooting.md](references/reporting-and-troubleshooting.md)
  for method reconciliation, cohort gotchas, reviewer questions, common errors,
  installation notes, references, and related Skills.
- Run or adapt
  [ancestral_reconstruction.py](scripts/ancestral_reconstruction.py) for the
  provider's complete Python demonstration. The method-specific extracted
  examples are `codeml_asr.py`, `iqtree_ancestral.sh`, `iqtree_state.py`,
  `stochastic_mapping.R`, `continuous_trait_asr.R`, and `grasp_asr.sh`.
- Use [usage-guide.md](usage-guide.md) for prompt examples and expected agent
  behavior.
- Read [provenance.md](references/provenance.md) for origin, licensing, and
  normalization boundaries.

## Required Deliverables

- Validated input summary and exact node or lineage target.
- Method-selection rationale, model or prior specification, root treatment,
  software versions, and reproducible command or script.
- Reconstructed states with uncertainty at every reported node or site.
- Model-fit, convergence, and sensitivity evidence appropriate to the method.
- Explicit list of ambiguous and model-sensitive states.
- Reconciliation of materially different reasonable methods.
- Limitations on biological interpretation and, for resurrection, an
  alternative-construct plan rather than a single unqualified sequence.
