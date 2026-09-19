# Trajectory Inference - Usage Guide

## Overview

Trajectory inference orders cells along a continuous differentiation manifold and assigns a pseudotime, fate probabilities, or RNA velocity direction. The central caution: a snapshot is not a movie. Pseudotime is geometry, not a clock; the existence of a continuum is a biological judgment made before ordering; the root choice flips every gene trend; and near bifurcations the honest output is a distribution over fates, not a hard branch label.

## Prerequisites

See SKILL.md's "Installation" section for R/Python install commands and the Windows/Monocle3 caveat.

## Quick Start

Tell your AI agent what you want to do:
- "Check whether my cells form a real trajectory or are discrete types"
- "Order cells by pseudotime from the stem-cell population"
- "Compute fate probabilities at the branch point"
- "Run RNA velocity and tell me whether it is trustworthy here"

## Example Prompts

### Topology and Continuum Test
> "Run PAGA and tell me whether these clusters are actually connected"
> "Is this a real continuum or a mixture of discrete cell types?"
> "Build a PAGA-initialized UMAP so the global topology is faithful"

### Pseudotime and Rooting
> "Order cells by diffusion pseudotime rooted at the HSC marker, not by eye"
> "Show how the gene trends change if the root cluster changes"
> "Which method fits a bifurcating tree best here?"

### Fate Probabilities
> "Give me fate probabilities at the branch point, not hard branch labels"
> "Run CellRank 2 with a pseudotime kernel and find terminal states"
> "Compute differentiation potential (entropy) across the manifold"

### RNA Velocity
> "Run scVelo dynamical mode and check velocity confidence first"
> "The velocity arrows point backward in my mature cells - is velocity valid here?"
> "Re-run velocity with a second quantifier and check the direction is stable"

## What the Agent Will Do

See SKILL.md's Governing Principle (5 rules) and Method Decision Table -- they define the decision sequence (continuum check, topology-first method choice, root anchoring, fate probabilities near bifurcations, velocity/quantifier validation) this Skill follows.

## Related Skills

single-cell/clustering - Leiden clusters and the kNN graph PAGA, DPT, and velocity moments depend on
single-cell/preprocessing - normalization, HVG, and PCA choices the inferred axis inherits
single-cell/lineage-tracing - orthogonal lineage ground truth that tests state-based fate calls
single-cell/cell-communication - downstream signaling analysis along the trajectory
differential-expression/deseq2-basics - pseudobulk DE between trajectory endpoints or branches
