# Structural Alignment - Usage Guide

## Overview

This skill covers backbone-aware structural alignment: pairwise (TM-align, US-align, Bio.PDB.Superimposer), database search at scale (Foldseek 3Di, DALI), and structural multiple alignment (Foldmason, MUSTANG). Use it when sequence identity drops below ~25% (twilight zone), when remote-homology detection is required, or when fold-similarity quantification is the goal. Install commands, thresholds and pitfalls are in `SKILL.md` (Install, metric tables, Common Errors).

## Example Prompts

### Pairwise Structural Alignment
> "Align these two crystal structures with TM-align and report the TM-score and superposition"

> "I have two homologs at 18% sequence identity. Compute structural alignment with US-align and tell me whether they share a fold."

> "Compare apo and holo conformations of the same protein and quantify the conformational change."

> "Score my AlphaFold model against the crystal structure: TM-score, RMSD and GDT-TS."

### Database Search at Scale
> "Search my predicted structure against AlphaFoldDB to find evolutionary distant homologs"

> "I have a metagenomic protein with no PDB hit. Predict its structure with ESMFold and search the ESM Atlas."

> "Cluster all PDB structures sharing TM-score > 0.5 with my query."

> "Search this antibody-antigen complex against AFDB-Multimer with Foldseek-Multimer and rank by complex TM-score"

### Structural Multiple Alignment
> "Build an MSA from these 30 PDB structures using Foldmason"

> "I need a structural MSA for tree inference from a remote homologue family"

### Twilight-Fold Homology
> "These two folds have TM-score 0.4. Re-align them with DALI and give me the Z-score."

## Related Skills

- alignment/multiple-alignment - Sequence MSA tools (MAFFT, MUSCLE5, ClustalOmega) for >25% identity
- alignment/pairwise-alignment - Sequence pairwise alignment when twilight zone is not crossed
- alignment/msa-statistics - Per-column statistics on structurally-derived MSAs
- alignment/alignment-trimming - Trim structural MSAs before phylogenetics
- phylogenetics/modern-tree-inference - Build trees from structural MSAs
