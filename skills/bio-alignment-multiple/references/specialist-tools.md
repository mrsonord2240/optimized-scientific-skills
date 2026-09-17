# Specialist Alignment Tools: BAli-Phy, OMM_MACSE, HyPhy, vcMSA

Moved out of `SKILL.md` (bio-alignment-multiple) because these apply only to joint alignment+tree
estimation, codon-pipeline wrapping for selection-grade input, or protein-language-model alignment at
low identity -- consult this file when the Critical Concepts, Codon-Aware Alignment, or Beyond MAFFT
and MUSCLE sections in `SKILL.md` point you here.

## Joint MSA-Phylogeny Co-estimation (Small Datasets Only)

The theoretically correct answer to guide-tree dependency is to estimate the alignment and tree jointly under a statistical evolutionary model rather than treating MSA as a fixed input to phylogenetic inference. BAli-Phy version 3 (Redelings 2021 Bioinf 37:3032) does this via MCMC, producing a posterior distribution over alignments and trees with insertion/deletion rates as model parameters. Version 3 is O(n) instead of O(n^2) per likelihood evaluation, but the practical ceiling remains ~70-200 sequences before runtime becomes prohibitive (weeks for >100 sequences). When the dataset fits the cap, BAli-Phy gives the most defensible alignment+tree pair for publication; when it does not, the practical alternative is MUSCLE5 ensemble + IQ-TREE per-replicate (Edgar 2022) to approximate posterior alignment uncertainty without joint estimation.

| Dataset size | Recommended approach |
|--------------|----------------------|
| < 70 sequences | BAli-Phy v3 joint MSA+tree posterior; gold standard |
| 70 - 200 sequences | BAli-Phy v3 if compute allows (weeks); else MUSCLE5 ensemble + IQ-TREE per replicate |
| > 200 sequences | MUSCLE5 ensemble (`-stratified`) + IQ-TREE per replicate; BAli-Phy not feasible |
| > 1000 sequences | Single MAFFT/MUSCLE5 + standard bootstrap; ensemble methods become intractable |

Redelings BD. 2021. BAli-Phy version 3: model-based co-estimation of alignment and phylogeny. Bioinf 37:3032-3034.

## OMM_MACSE: Recommended Pipeline Wrapper

OMM_MACSE (Ranwez group, used in OrthoMaM v10+) chains: MACSE `trimNonHomologousFragments` to remove non-homologous fragments (annotation errors, retained introns/UTRs), MAFFT pre-alignment for guide-tree, MACSE v2 frameshift-aware refinement, then soft HMMcleaner cleaning. This is the recommended pipeline for genome-scale ortholog datasets where some genes will contain frameshifts.

```bash
OMM_MACSE_v12.02.sif --in_seq_file orthogroup.fasta --out_dir omm_out --out_file_prefix orthogroup --genetic_code_number 1
```

## HyPhy pre-msa.bf / post-msa.bf

For HyPhy-grade dN/dS analyses (BUSTED, MEME, aBSREL, RELAX), Pond lab's standard workflow is:

```bash
hyphy pre-msa.bf --input cds.fasta
mafft --auto cds.fasta_protein.fas > cds.fasta_protein.msa
hyphy post-msa.bf --protein-msa cds.fasta_protein.msa --nucleotide-sequences cds.fasta_nuc.fas --output cds.codon.msa
```

`pre-msa.bf` strips internal stop codons and translates; `post-msa.bf` validates that all sequences remain in-frame after threading. Without this validation step, a single mis-aligned codon can cascade into a spurious episodic-selection call.

## vcMSA: Mixed-Length Input Limitation

**vcMSA limitation:** Pre-filter input to within ~2x mean length; vcMSA degrades sharply on mixed full-length/fragment input because ProtT5 embeddings encode positional context. For mixed-length sets, segment long sequences via HHsearch domain decomposition before alignment, or use Foldmason on predicted structures.
