---
name: bio-single-cell-splicing
category: Data Analysis
description: >-
  Analyzes alternative splicing at single-cell resolution. Start by determining
  whether library chemistry supports the requested inference: 10X 3' is
  unsuitable for transcriptome-wide splicing; full-length plate and single-cell
  long-read assays support per-cell isoform structure. Routes work to MARVEL,
  BRIE2, scQuint, SpliZ, Psix, Sierra (APA), or replicate-aware pseudobulk.
tool_type: mixed
primary_tool: MARVEL
license: MIT
author: GPTomics
---

# Single-Cell Splicing Analysis

The first decision is **chemistry**, not tool. Most droplet 3' scRNA-seq cannot support transcriptome-wide splicing inference: reverse transcription starts at poly(A), leaving most informative internal junctions unobserved. Plate-based full-length methods and single-cell long-read sequencing recover gene-body isoform structure.

## Version Compatibility

Reference examples were checked with MARVEL 2.0.5 (R 4.4, Seurat 5.5), BRIE2 2.3.0, scQuint 0.3.3, Psix 0.10.8, Sierra 0.99.27, regtools 1.0.0, leafcutter 0.2.9, anndata 0.12, scanpy 1.11, and pandas 2.3. SpliZ was not run because Nextflow was unavailable. Installation commands in the tool references are guidance, not a claim that installation was re-run in your environment.

Before copying a pattern, verify versions and the local API:

```bash
pip show <package>
<tool> --version
<tool> --help
```

```r
packageVersion('<pkg>')
?function_name
```

If code raises `ImportError`, `AttributeError`, or `TypeError`, inspect the installed signature and adapt it; do not retry an incompatible call unchanged.

## Chemistry Gate

| Chemistry | Transcriptome-wide splicing? | Appropriate route |
|---|---|---|
| 10X 3' (Chromium v3, GEM-X, v4, Flex) | No; a few near-3' events may be observable | Sierra for APA, or generate full-length/long-read data |
| 10X 5' GEX | Limited to 5'-proximal events | Alternative TSS/APA question, or MAS-Iso-seq |
| Smart-seq2 | Yes | MARVEL or BRIE2 |
| Smart-seq3 / Smart-seq3xpress / FLASH-seq | Yes, with UMI-aware molecule counts where available | MARVEL or BRIE2 |
| VASA-seq / STORM-seq | Yes; total RNA can expose intron retention | MARVEL, interpreting IR carefully |
| MAS-Iso-seq + 10X 5', scISOr-Seq2, ONT scRNA | Yes; full isoforms | FLAMES, IsoQuant, or long-read-splicing |

10X 5' GEX is not a transcriptome-wide splicing solution: it shifts capture toward the 5' end. V(D)J recovery requires the dedicated 10X Immune Profiling kit with TCR/BCR enrichment, not 5' GEX alone.

For 10X 3', do not promise per-cell cassette-exon PSI or cell-type-specific transcriptome-wide exon calls. Explain the chemistry limitation, determine whether the actual question is APA, and offer Sierra, replicate-aware pseudobulk only where junction support exists, or a full-length/long-read design.

## Choose a Route

| Goal | Route | Detailed workflow |
|---|---|---|
| Smart-seq cassette exons and cell-type comparisons | MARVEL | [MARVEL](references/marvel.md) |
| Per-cell PSI with uncertainty and covariate testing | BRIE2 | [BRIE2](references/brie2.md) |
| Annotation-free, plate-based junction clusters | scQuint | [scQuint](references/scquint.md) |
| Annotation-free cell-state discovery | SpliZ | [SpliZ](references/spliz.md) |
| Regulated AS over a trajectory | Psix | [Psix](references/psix.md) |
| 10X 3' alternative polyadenylation | Sierra | [Sierra](references/sierra.md) |
| Differential splicing between cell types with donors/plates | Pseudobulk leafcutter or rMATS | [Pseudobulk](references/pseudobulk.md) |
| Full-length isoforms per cell | MAS-Iso-seq + FLAMES / IsoQuant | long-read-splicing |

scQuint is validated for plate-based data and its authors recommend against 10X 3'/5' analysis because capture bias confounds junction usage. Sierra is APA, not cassette-exon splicing. SpliZ is an annotation-free score and can reveal a gene-level splicing signal without yielding a conventional PSI event.

## Hard Boundaries

- Do not impute PSI with expression imputation tools (MAGIC, scImpute, ALRA): averaging missing junction evidence erases the heterogeneity being tested. Psix tests smoothness of observed PSI instead.
- Do not treat cells as independent biological replicates for between-condition inference. Sum junctions by **sample x cell type** and require at least three donor, batch, or plate replicates per group.
- Do not interpret intron retention in snRNA-seq as mature-isoform regulation without accounting for nuclear, incompletely spliced RNA.
- Filter doublets before splicing analysis; mixed cells can generate artificial middle or bimodal PSI distributions.
- Use the same genome build and chromosome naming convention across alignments, event annotation, and downstream tables.

## Interpretation and QC

| Metric | Recommendation |
|---|---|
| Cells supporting an event | At least 50 for per-cell conclusions |
| Junction reads per event per cell | 5-10 for stable PSI; <=1 is unreliable |
| Per-cell PSI variance | Look for low within-cluster and large between-cluster separation, but inspect coverage first |
| Pseudobulk replicates | >=3 samples per group; >=50 cells per cluster is a practical minimum |
| Library choice | Full-length plate or long-read for transcriptome-wide AS; 3' chemistry for APA |

MARVEL modality labels are descriptive, not proof of mechanism: included (PSI near 1), excluded (near 0), bimodal (mixed states or bursting), middle (often technical mixture/low coverage), and multimodal. Confirm unusual modality with full-length data and doublet checks.

Beta-binomial or Dirichlet-multinomial models help with sparse overdispersed counts; they do not make unsupported 3' per-cell PSI reliable. Pseudobulk trades within-cluster heterogeneity for statistical power.

## Reconcile Disagreement

| Pattern | Likely cause | Next action |
|---|---|---|
| MARVEL significant, BRIE2 not | BRIE2 is more conservative about per-cell uncertainty | Use MARVEL for cell-type comparisons; retain BRIE2 for within-cluster uncertainty |
| BRIE2 significant, MARVEL not | Smooth cell-state association, not a discrete boundary | Test the trajectory with Psix |
| SpliZ significant, MARVEL not | Novel or differently represented junction structure | Inspect junctions manually |
| Sierra significant, MARVEL not | APA and splicing are different biology | Report them separately |
| Pseudobulk significant, per-cell method not | Aggregate power exceeds per-cell power | Report a cluster-level result |

## Failure-Mode Index

| Symptom | Cause | Action |
|---|---|---|
| MARVEL `$ operator is invalid for atomic vectors` or all-NA PSI | Long rather than wide SJ matrix, missing `coord.intron`, contig mismatch, or too-high coverage threshold | Follow [MARVEL input checks](references/marvel.md#input-contract-and-ri) |
| MARVEL RI errors with `argument is of length zero` | RI needs `IntronCounts`, `thread >= 2`, and an explicit read length | Use the RI call in [MARVEL](references/marvel.md#input-contract-and-ri) |
| `brie-count` cannot retrieve an index | BAM lacks `.bai` | `samtools index` every BAM |
| BRIE2 output has unexpectedly few events | Silent count/cell/MIF filters | Inspect the log and lower filters only with a justified coverage threshold |
| scQuint produces zero annotated introns | scQuint prepends `chr` to GTF contigs | Apply the explicit [scQuint contig contract](references/scquint.md#contig-contract) and check `adata.n_vars` |
| scQuint sparse groups error instead of returning tables | All groups were filtered | Catch the `ValueError`, lower both thresholds only for adequately covered plate data |
| Psix constructor/query error | It expects PSI and mRNA tables plus a latent space; results use `qvals` | Use [Psix](references/psix.md) |
| Sierra `FindPeaks` says unused argument | Parameter is `bamfile`, and `junctions.file` is required | Use [Sierra](references/sierra.md) |

## Related Skills

- single-cell/preprocessing - QC and normalization
- single-cell/clustering - Cell-type annotation prerequisite
- single-cell/doublet-detection - doublet filtering
- single-cell/data-io - h5ad / Seurat I/O
- splicing-quantification - bulk RNA-seq context
- long-read-splicing - full-isoform analysis from MAS-Iso-seq and scISOr-Seq2

## References

- Huang & Sanguinetti 2021, *Genome Biology* - BRIE2
- Wen et al. 2023, *Nucleic Acids Research* 51:e29 - MARVEL
- Benegas, Fischer & Song 2022, *eLife* - scQuint
- Olivieri et al. 2022, *Nature Methods* - SpliZ
- Buen Abad Najar et al. 2022, *Genome Research* 32:1385 - Psix
- Patrick et al. 2020, *Genome Biology* - Sierra
- Song et al. 2017, *Molecular Cell* - splicing modality
- Hagemann-Jensen et al. 2020/2022, *Nature Biotechnology* - Smart-seq3/3xpress
- Hahaut et al. 2022, *Nature Biotechnology* - FLASH-seq
- Salmen et al. 2022, *Nature Biotechnology* - VASA-seq
- Al'Khafaji et al. 2024, *Nature Biotechnology* - MAS-Iso-seq / Kinnex
- Joglekar et al. 2024, *Nature Neuroscience* 27:1051-1063 - scISOr-Seq2
