# Cell-Cell Communication - Usage Guide

## Overview

Cell-cell communication (CCC) analysis ranks ligand-receptor interactions between cell types from scRNA-seq. Every output is a co-expression proxy, not proof of signaling, and competing methods disagree because they estimate different quantities (specificity vs magnitude vs probability). The defensible workflow is consensus-first (LIANA), with a resource-sensitivity check and orthogonal validation, reserving NicheNet for the distinct mechanistic question of which ligand drives a receiver's response. See SKILL.md's Prerequisites section for install commands.

## Quick Start

Tell your AI agent what you want to do:
- "Rank ligand-receptor interactions between my cell types with a consensus method"
- "Check whether my top interactions survive a different L-R database"
- "Which ligand from macrophages best explains the DE genes in T cells?"
- "Compare cell communication between control and disease"
- "Summarize WNT signaling at the pathway level with sender and receiver roles"

## Example Prompts

### Consensus Inference
> "Run LIANA rank_aggregate and report both magnitude and specificity ranks"
> "Find robust ligand-receptor pairs between fibroblasts and epithelial cells"

### Resource Sensitivity
> "Re-run the same method with the CellPhoneDB and CellChat resources and show which pairs are stable"
> "Is my VEGF-VEGFR hit robust to the choice of database?"

### Specificity Testing
> "Run CellPhoneDB with permutation p-values and proper complex handling"
> "Use the DEG-based CellPhoneDB method instead of one-vs-rest"

### Pathway-Level and Roles
> "Summarize signaling pathways between these cell types with CellChat"
> "Which cell types are the dominant senders and receivers of TGF-beta?"

### Downstream Mechanism (NicheNet)
> "Which ligands explain the activated-T-cell gene signature?"
> "Rank ligands by downstream target-gene activity in the receiver"

### Condition Comparison
> "Compare communication between control and treatment and show gained/lost interactions"

## Related Skills

single-cell/cell-annotation - Cell-type labels define senders and receivers; annotation resolution is a hidden CCC hyperparameter
single-cell/clustering - Cluster granularity changes who is "specific"; fix it before running CCC
single-cell/doublet-detection - Doublets create fake co-expressing cells that masquerade as senders-receivers
single-cell/preprocessing - Ambient-RNA decontamination and stress-gene handling happen here, before CCC
single-cell/metabolite-communication - Metabolite-mediated CCC (enzyme-sensor) as the doubly-inferred counterpart to ligand-receptor
spatial-transcriptomics/spatial-communication - Proximity-constrained CCC when spatial coordinates are available
pathway-analysis/go-enrichment - Functional enrichment of NicheNet target genes or interacting receptors
differential-expression/deseq2-basics - Pseudobulk DE to build the receiver gene set NicheNet requires
