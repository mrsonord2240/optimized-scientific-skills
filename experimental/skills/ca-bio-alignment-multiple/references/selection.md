# Selection



When a tool is failing on a dataset, switch to a tool from a different algorithmic family rather than tuning flags. The six families and their characteristic failure modes:

| Family | Representative tools | Best at | Fails when |
|--------|---------------------|---------|-----------|
| Progressive | ClustalW, MAFFT FFT-NS-2 | Fast, large datasets, similar lengths | Early-stage gap errors propagate; no recovery |
| Iterative refinement | MAFFT L-INS-i, MUSCLE3, PRRN | Recovers from progressive errors at <2000 seqs | Slow on >2000; still guide-tree dependent |
| Consistency-based | T-Coffee, ProbCons | Highest accuracy <100 seqs; integrates evidence | O(N^2 to N^4) scaling; heavy compute |
| HMM-based | HMMER hmmalign, ClustalOmega (HHalign), UPP, WITCH | Adding sequences to a curated profile; fragmentary input | Needs an existing high-quality profile or backbone |
| Divide-and-conquer | PASTA, MAGUS, MUSCLE5 super5 | Heterogeneous large datasets (>10k seqs) | Sub-alignment merges can introduce artefacts |
| Structure or pLM-informed | Foldmason, PROMALS3D, vcMSA, 3D-Coffee | Dark proteome, <15% identity, dataset has structures | Requires structures or a working pLM |

## Tool Selection

Pick by dataset size and divergence; the default recommendations follow the table.

| Tool | Best For | Max Sequences | Accuracy | Speed |
|------|----------|---------------|----------|-------|
| MAFFT L-INS-i | Highest accuracy, <200 seqs | ~200 | Highest | Slow |
| MAFFT FFT-NS-2 | Large datasets, good balance | ~50,000 | Good | Fast |
| MAFFT E-INS-i | Sequences with long unalignable internal regions | ~200 | High | Slow |
| MUSCLE5 (PPP `-align`) | Benchmarked highest accuracy on Balifam-10000 | ~1000 | Highest | Medium |
| MUSCLE5 (`-super5`) | Large datasets via mBed clustering | ~100,000+ | Good | Medium |
| ClustalOmega | Very large datasets, HMM-based profiles | ~190,000 in published benchmark (Sievers et al 2011 Mol Syst Biol) | Good | Fast |
| T-Coffee (default) | Small datasets needing maximum accuracy | ~200 | Highest | Slowest |

**Default recommendation**: MAFFT L-INS-i for <200 sequences; MAFFT FFT-NS-2 or MUSCLE5 super5 for thousands; MUSCLE5 ensemble (`-stratified`) when alignment confidence estimates are needed.

### Cross-Aligner Sensitivity Check

When downstream analysis depends on a specific column (a candidate selection site, a contact-prediction position), run BOTH MAFFT L-INS-i and MUSCLE5 `-align` and verify that column is stable across the two outputs. If unstable, flag as low-confidence regardless of GUIDANCE2/TCS scores. The MUSCLE5 ensemble (`-stratified`/`-diversified`) accomplishes the same check directly within one tool and is preferred when ensemble output is acceptable downstream.

### Beyond MAFFT and MUSCLE: Scale and Domain

Some workloads exceed what the four main tools handle gracefully. Use the scale-and-domain table below to escape default-tool blind spots.

| Scenario | Recommended tool | Why |
|----------|------------------|-----|
| 100k - millions of sequences (UniRef cluster reps) | FAMSA v2 (Deorowicz et al 2016 SREP, v2 2024) | Million-scale MSA in hours with Pareto-optimal accuracy |
| Heterogeneous large dataset (variable length, divergence) | PASTA (Mirarab et al 2015 J Comp Biol) or MAGUS (Smirnov & Warnow 2021 Bioinf) | Divide-and-conquer plus iterative re-alignment |
| Fragmentary input (metagenomics, eDNA, ancient DNA) | UPP (Nguyen et al 2015 Genome Biol) or WITCH (Shen et al 2022 J Comp Biol) | HMM backbone tolerates partial sequences |
| Dark proteome / <15% identity | vcMSA (McWhite, Armour-Garb & Singh 2023 Genome Res; alpha-stage tool per its repo README; last release Oct 2023; test on representative inputs before pipeline use) or Foldmason (Gilchrist et al 2024) | pLM-embedding or structural alignment where sequence fails |
| RNA with secondary structure | Infernal `cmalign` (Nawrocki & Eddy 2013), R-Coffee, MAFFT-Q-INS-i | Sequence + base-pair consistency |
| Strand-unknown (Sanger, Nanopore raw) | `mafft --adjustdirection --globalpair` | Auto-detects and reverse-complements as needed |
| Adding many sequences to a curated profile | HMMER `hmmalign` (Eddy 2011 PLOS CB) | Profile-driven; better than `mafft --add` for Pfam-style use |

**vcMSA limitation:** mixed full-length/fragment input degrades vcMSA sharply because ProtT5 embeddings encode positional context. See `references/specialist-tools.md` for the pre-filtering and segmentation workaround.

**`--adjustdirection` renames reversed sequences.** MAFFT silently prepends `_R_` to the header of every sequence it reverse-complements (verified on MAFFT 7.526: `mafft --adjustdirection --globalpair small.fa` turned `>revcomp_seq` into `>_R_revcomp_seq`, leaving already-forward sequences untouched). Strip the `_R_` prefix from sequence IDs before any downstream join or merge keyed on sequence name (tree tip labels, metadata tables, `pal2nal` pairing) -- an unstripped prefix silently breaks the key match instead of erroring.

## Critical Concepts

### Mitigate Guide-Tree Dependency

All major MSA tools build a guide tree, then align progressively along it; once a gap is inserted in the progressive phase, it is never removed, so early errors propagate. To mitigate: prefer iterative-refinement modes (MAFFT `-i`, MUSCLE5), use consistency scoring (T-Coffee) for small datasets, and quantify uncertainty with GUIDANCE2 bootstrapping or the MUSCLE5 ensemble before publishing column-specific conclusions.

### Joint MSA-Phylogeny Co-estimation (Small Datasets Only)

The theoretically correct answer to guide-tree dependency is to estimate the alignment and tree jointly under a statistical evolutionary model rather than treating MSA as a fixed input to phylogenetic inference, via BAli-Phy's MCMC co-estimation (practical ceiling ~70-200 sequences); the practical alternative above that ceiling is MUSCLE5 ensemble + IQ-TREE per-replicate (Edgar 2022). See `references/specialist-tools.md` for the BAli-Phy workflow and the dataset-size table.

### Sequence Divergence Thresholds

| Protein Identity | Signal Level | Recommendation |
|-----------------|--------------|----------------|
| >40% | Strong | Any MSA tool produces reliable alignment |
| 25-40% | Moderate (twilight zone begins) | Use iterative methods (L-INS-i, MUSCLE5); validate with GUIDANCE2 |
| 20-25% | Weak | Profile-profile methods (HHpred); consider structural alignment |
| <15-20% (length-dependent twilight) | Noise dominates signal | Sequence MSA is unreliable; switch to structural alignment (Foldseek, TM-align) or pLM aligners -- see `bio-alignment-structural` |
