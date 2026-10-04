# Remaining Skills

Not yet refined. Scope is deliberately limited to the rest of
[GPTomics/bioSkills](https://github.com/GPTomics/bioSkills) at commit
`d91ed3d563019e649dc854c56ccd62551359488a`; other source corpora are out of scope for now.

**366 remaining** across 48 folders. 195 are already refined and live in `skills/`.

The source tree holds 562 Skills: 195 refined, 0 audited and excluded, 1 out of scope, 366 remaining.

Folders are ordered by `usage`, the mean score of their remaining Skills. A score is
0-100: half conda downloads of the Skill's primary tool, half how many workflow Skills
depend on it. It is a proxy from `tools/usage_rank.py` in the records repository, not
measured use.

| folder | remaining | refined | usage |
| --- | ---: | ---: | ---: |
| read-alignment | 4 | 0 | 56 |
| rna-quantification | 4 | 0 | 51 |
| phasing-imputation | 4 | 0 | 49 |
| read-qc | 7 | 0 | 49 |
| variant-calling | 7 | 6 | 49 |
| genome-assembly | 9 | 0 | 48 |
| differential-expression | 5 | 1 | 48 |
| liquid-biopsy | 5 | 2 | 47 |
| hi-c-analysis | 9 | 0 | 47 |
| multi-omics-integration | 5 | 0 | 47 |
| epidemiological-genomics | 5 | 0 | 46 |
| small-rna-seq | 6 | 0 | 46 |
| clinical-biostatistics | 11 | 1 | 45 |
| flow-cytometry | 8 | 0 | 45 |
| population-genetics | 6 | 1 | 45 |
| metagenomics | 8 | 0 | 45 |
| genome-annotation | 7 | 0 | 44 |
| ribo-seq | 6 | 0 | 44 |
| restriction-analysis | 5 | 0 | 43 |
| sequence-io | 8 | 1 | 43 |
| sequence-manipulation | 6 | 1 | 43 |
| expression-matrix | 5 | 0 | 41 |
| epitranscriptomics | 5 | 0 | 41 |
| imaging-mass-cytometry | 7 | 0 | 40 |
| reporting | 6 | 0 | 39 |
| genome-intervals | 8 | 0 | 38 |
| structural-biology | 10 | 0 | 38 |
| systems-biology | 7 | 0 | 38 |
| workflow-management | 5 | 0 | 36 |
| genome-engineering | 5 | 0 | 35 |
| long-read-sequencing | 9 | 0 | 35 |
| temporal-genomics | 6 | 0 | 35 |
| methylation-analysis | 10 | 0 | 35 |
| ecological-genomics | 6 | 0 | 35 |
| chip-seq | 11 | 1 | 34 |
| data-visualization | 10 | 10 | 34 |
| spatial-transcriptomics | 12 | 0 | 34 |
| copy-number | 11 | 0 | 33 |
| rna-structure | 4 | 0 | 33 |
| clinical-databases | 7 | 5 | 30 |
| workflows | 37 | 4 | 29 |
| tcr-bcr-analysis | 6 | 0 | 26 |
| primer-design | 4 | 0 | 25 |
| immunoinformatics | 6 | 0 | 24 |
| gene-regulatory-networks | 6 | 0 | 24 |
| clip-seq | 11 | 1 | 22 |
| comparative-genomics | 12 | 1 | 20 |
| chemoinformatics | 5 | 15 | 15 |

## Promoted, fix pass still needed

These Skills are in `skills/` because their audit found them deployable with no open P0.
They have not yet been through a fix pass, however high they scored. Their open findings
are in the audit record. `fix_pass` in PROVENANCE.json carries the same flag.

| skill | score | grade |
| --- | ---: | --- |

## Promoted, re-audit still needed

Changed in the provider after its latest audit, so the score describes earlier bytes.
`reaudit` in PROVENANCE.json carries the same flag.

| skill | score at last audit | audited on |
| --- | ---: | --- |

## Audited and excluded

Audited and did not pass. Not pending — rejected until the defects behind the score are
fixed.

| skill | score | grade | open P0 |
| --- | ---: | --- | ---: |

## Out of scope

| skill | reason |
| --- | --- |
| `clawhub-installer` | upstream's own corpus installer, not a science Skill; declares os: darwin/linux only and exists to install the other Skills |

## The list

### chemoinformatics

- `bio-reaction-enumeration` — `chemoinformatics/reaction-enumeration` — usage 43 (RDKit)
- `bio-free-energy-calculations` — `chemoinformatics/free-energy-calculations` — usage 34 (OpenFE)
- `bio-generative-design` — `chemoinformatics/generative-design` — usage 0 (REINVENT)
- `bio-ml-docking-rescoring` — `chemoinformatics/ml-docking-rescoring` — usage 0 (DiffDock)
- `bio-retrosynthesis` — `chemoinformatics/retrosynthesis` — usage 0 (AiZynthFinder)

### chip-seq

- `bio-chipseq-qc` — `chip-seq/chipseq-qc` — usage 53 (deepTools, in 1 workflow)
- `bio-chipseq-visualization` — `chip-seq/chipseq-visualization` — usage 53 (deepTools, in 1 workflow)
- `bio-chipseq-differential-binding` — `chip-seq/differential-binding` — usage 48 (DiffBind, in 1 workflow)
- `bio-chipseq-motif-analysis` — `chip-seq/motif-analysis` — usage 47 (HOMER, in 1 workflow)
- `bio-chipseq-peak-annotation` — `chip-seq/peak-annotation` — usage 46 (ChIPseeker, in 1 workflow)
- `bio-chipseq-peak-calling` — `chip-seq/peak-calling` — usage 44 (macs3, in 1 workflow)
- `bio-chipseq-spike-in-normalization` — `chip-seq/spike-in-normalization` — usage 34 (DiffBind)
- `bio-chipseq-chromatin-state-segmentation` — `chip-seq/chromatin-state-segmentation` — usage 29 (ChromHMM)
- `bio-chipseq-cut-and-run-tag` — `chip-seq/cut-and-run-tag` — usage 25 (SEACR)
- `bio-chipseq-chip-deep-learning` — `chip-seq/chip-deep-learning` — usage 0 (chrombpnet)
- `bio-chipseq-super-enhancers` — `chip-seq/super-enhancers` — usage 0 (ROSE)

### clinical-biostatistics

- `bio-clinical-biostatistics-categorical-tests` — `clinical-biostatistics/categorical-tests` — usage 65 (scipy, in 1 workflow)
- `bio-clinical-biostatistics-effect-measures` — `clinical-biostatistics/effect-measures` — usage 62 (statsmodels, in 1 workflow)
- `bio-clinical-biostatistics-logistic-regression` — `clinical-biostatistics/logistic-regression` — usage 62 (statsmodels, in 1 workflow)
- `bio-clinical-biostatistics-power-sample-size` — `clinical-biostatistics/power-and-sample-size` — usage 62 (statsmodels, in 1 workflow)
- `bio-clinical-biostatistics-subgroup-analysis` — `clinical-biostatistics/subgroup-analysis` — usage 62 (statsmodels, in 1 workflow)
- `bio-clinical-biostatistics-cdisc-data` — `clinical-biostatistics/cdisc-data-handling` — usage 54 (pyreadstat, in 1 workflow)
- `bio-clinical-biostatistics-survival-analysis` — `clinical-biostatistics/survival-analysis` — usage 52 (lifelines, in 1 workflow)
- `bio-clinical-biostatistics-trial-reporting` — `clinical-biostatistics/trial-reporting` — usage 49 (tableone, in 1 workflow)
- `bio-clinical-biostatistics-missing-data` — `clinical-biostatistics/missing-data-sensitivity` — usage 15 (rbmi, in 1 workflow)
- `bio-clinical-biostatistics-multiplicity-graphical` — `clinical-biostatistics/multiplicity-graphical` — usage 15 (gMCP, in 1 workflow)
- `bio-clinical-biostatistics-bayesian-trials` — `clinical-biostatistics/bayesian-trials` — usage 0 (RBesT)

### clinical-databases

- `bio-clinical-databases-variant-prioritization` — `clinical-databases/variant-prioritization` — usage 50 (pandas)
- `bio-clinical-databases-hla-typing` — `clinical-databases/hla-typing` — usage 42 (T1K, in 1 workflow)
- `bio-clinical-databases-tumor-mutational-burden` — `clinical-databases/tumor-mutational-burden` — usage 38 (cyvcf2)
- `bio-clinical-databases-msi-detection` — `clinical-databases/msi-detection` — usage 32 (MSIsensor-pro)
- `bio-clinical-databases-somatic-signatures` — `clinical-databases/somatic-signatures` — usage 24 (SigProfilerAssignment)
- `bio-clinical-databases-polygenic-risk` — `clinical-databases/polygenic-risk` — usage 23 (PGS Catalog Calculator)
- `bio-clinical-databases-pharmacogenomics` — `clinical-databases/pharmacogenomics` — usage 0 (PharmCAT)

### clip-seq

- `bio-clip-seq-clip-alignment` — `clip-seq/clip-alignment` — usage 39 (STAR)
- `bio-clip-seq-clip-preprocessing` — `clip-seq/clip-preprocessing` — usage 34 (umi_tools)
- `bio-clip-seq-clip-motif-analysis` — `clip-seq/clip-motif-analysis` — usage 32 (HOMER)
- `bio-clip-seq-binding-site-annotation` — `clip-seq/binding-site-annotation` — usage 31 (ChIPseeker)
- `bio-clip-seq-clip-qc` — `clip-seq/clip-qc` — usage 30 (preseq)
- `bio-clip-seq-crosslink-site-detection` — `clip-seq/crosslink-site-detection` — usage 29 (PureCLIP)
- `bio-clip-seq-differential-clip` — `clip-seq/differential-clip` — usage 27 (DEWSeq)
- `bio-clip-seq-clip-peak-calling` — `clip-seq/clip-peak-calling` — usage 21 (CLIPper)
- `bio-clip-seq-clip-deep-learning` — `clip-seq/clip-deep-learning` — usage 0 (RBPNet)
- `bio-clip-seq-m6a-clip` — `clip-seq/m6a-clip` — usage 0 (miCLIP2)
- `bio-clip-seq-stamp-antibody-free` — `clip-seq/stamp-antibody-free` — usage 0 (STAMP)

### comparative-genomics

- `bio-comparative-genomics-positive-selection` — `comparative-genomics/positive-selection` — usage 35 (PAML)
- `bio-comparative-genomics-ortholog-inference` — `comparative-genomics/ortholog-inference` — usage 32 (OrthoFinder)
- `bio-comparative-genomics-genome-distance-and-species-delineation` — `comparative-genomics/genome-distance-and-species-delineation` — usage 31 (skani)
- `bio-comparative-genomics-pangenome-analysis` — `comparative-genomics/pangenome-analysis` — usage 30 (Panaroo)
- `bio-comparative-genomics-gene-family-evolution` — `comparative-genomics/gene-family-evolution` — usage 28 (CAFE5)
- `bio-comparative-genomics-whole-genome-alignment` — `comparative-genomics/whole-genome-alignment` — usage 28 (Cactus)
- `bio-comparative-genomics-whole-genome-duplication` — `comparative-genomics/whole-genome-duplication` — usage 28 (wgd)
- `bio-comparative-genomics-synteny-analysis` — `comparative-genomics/synteny-analysis` — usage 23 (MCScanX)
- `bio-comparative-genomics-comparative-annotation-projection` — `comparative-genomics/comparative-annotation-projection` — usage 0 (TOGA)
- `bio-comparative-genomics-gene-tree-species-tree-reconciliation` — `comparative-genomics/gene-tree-species-tree-reconciliation` — usage 0 (ALE)
- `bio-comparative-genomics-hgt-detection` — `comparative-genomics/hgt-detection` — usage 0 (HGTector)
- `bio-comparative-genomics-introgression-detection` — `comparative-genomics/introgression-detection` — usage 0 (Dsuite)

### copy-number

- `bio-copy-number-cnv-visualization` — `copy-number/cnv-visualization` — usage 63 (matplotlib, in 1 workflow)
- `bio-copy-number-cnvkit-analysis` — `copy-number/cnvkit-analysis` — usage 58 (cnvkit, in 2 workflows)
- `bio-copy-number-cnv-annotation` — `copy-number/cnv-annotation` — usage 56 (bedtools, in 1 workflow)
- `bio-copy-number-gatk-cnv` — `copy-number/gatk-cnv` — usage 53 (gatk, in 1 workflow)
- `bio-copy-number-copy-ratio-segmentation` — `copy-number/copy-ratio-segmentation` — usage 49 (DNAcopy, in 1 workflow)
- `bio-copy-number-allele-specific-copy-number` — `copy-number/allele-specific-copy-number` — usage 43 (ascat, in 1 workflow)
- `bio-copy-number-focal-amplification-ecdna` — `copy-number/focal-amplification-ecdna` — usage 27 (AmpliconArchitect)
- `bio-copy-number-recurrent-cnv` — `copy-number/recurrent-cnv` — usage 15 (gistic2, in 1 workflow)
- `bio-copy-number-germline-cnv-interpretation` — `copy-number/germline-cnv-interpretation` — usage 0 (ClassifyCNV)
- `bio-copy-number-hrd-scoring` — `copy-number/hrd-scoring` — usage 0 (scarHRD)
- `bio-copy-number-subclonal-copy-number` — `copy-number/subclonal-copy-number` — usage 0 (battenberg)

### data-visualization

- `bio-data-visualization-interactive-visualization` — `data-visualization/interactive-visualization` — usage 45 (plotly)
- `bio-data-visualization-distribution-plots` — `data-visualization/distribution-plots` — usage 42 (ggplot2)
- `bio-data-visualization-multipanel-figures` — `data-visualization/multipanel-figures` — usage 36 (patchwork)
- `bio-data-visualization-circos-plots` — `data-visualization/circos-plots` — usage 35 (circlize)
- `bio-data-visualization-heatmaps-clustering` — `data-visualization/heatmaps-clustering` — usage 33 (ComplexHeatmap)
- `bio-data-visualization-flow-and-transition-plots` — `data-visualization/flow-and-transition-plots` — usage 32 (ggalluvial)
- `bio-data-visualization-genome-tracks` — `data-visualization/genome-tracks` — usage 31 (pyGenomeTracks)
- `bio-data-visualization-lollipop-protein-maps` — `data-visualization/lollipop-protein-maps` — usage 30 (maftools)
- `bio-data-visualization-manhattan-qq-locuszoom` — `data-visualization/manhattan-qq-locuszoom` — usage 30 (qqman)
- `bio-data-visualization-upset-plots` — `data-visualization/upset-plots` — usage 30 (ComplexUpset)

### differential-expression

- `bio-differential-expression-de-results` — `differential-expression/de-results` — usage 52 (DESeq2, in 1 workflow)
- `bio-differential-expression-de-visualization` — `differential-expression/de-visualization` — usage 52 (DESeq2, in 1 workflow)
- `bio-differential-expression-timeseries-de` — `differential-expression/timeseries-de` — usage 52 (DESeq2, in 1 workflow)
- `bio-differential-expression-edger-basics` — `differential-expression/edger-basics` — usage 52 (edgeR, in 1 workflow)
- `bio-differential-expression-batch-correction` — `differential-expression/batch-correction` — usage 32 (sva)

### ecological-genomics

- `bio-ecological-genomics-community-ecology` — `ecological-genomics/community-ecology` — usage 54 (vegan, in 1 workflow)
- `bio-ecological-genomics-edna-metabarcoding` — `ecological-genomics/edna-metabarcoding` — usage 52 (dada2, in 1 workflow)
- `bio-ecological-genomics-biodiversity-metrics` — `ecological-genomics/biodiversity-metrics` — usage 44 (iNEXT, in 1 workflow)
- `bio-ecological-genomics-conservation-genetics` — `ecological-genomics/conservation-genetics` — usage 30 (hierfstat)
- `bio-ecological-genomics-landscape-genomics` — `ecological-genomics/landscape-genomics` — usage 29 (LEA)
- `bio-ecological-genomics-species-delimitation` — `ecological-genomics/species-delimitation` — usage 0 (ASAP)

### epidemiological-genomics

- `bio-epidemiological-genomics-amr-surveillance` — `epidemiological-genomics/amr-surveillance` — usage 48 (AMRFinderPlus, in 1 workflow)
- `bio-epidemiological-genomics-variant-surveillance` — `epidemiological-genomics/variant-surveillance` — usage 48 (Pangolin, in 1 workflow)
- `bio-epidemiological-genomics-pathogen-typing` — `epidemiological-genomics/pathogen-typing` — usage 47 (chewBBACA, in 1 workflow)
- `bio-epidemiological-genomics-transmission-inference` — `epidemiological-genomics/transmission-inference` — usage 45 (TransPhylo, in 1 workflow)
- `bio-epidemiological-genomics-phylodynamics` — `epidemiological-genomics/phylodynamics` — usage 44 (BEAST2, in 1 workflow)

### epitranscriptomics

- `bio-epitranscriptomics-merip-preprocessing` — `epitranscriptomics/merip-preprocessing` — usage 54 (STAR, in 1 workflow)
- `bio-epitranscriptomics-modification-visualization` — `epitranscriptomics/modification-visualization` — usage 43 (Guitar, in 1 workflow)
- `bio-epitranscriptomics-m6a-differential` — `epitranscriptomics/m6a-differential` — usage 41 (exomePeak2, in 1 workflow)
- `bio-epitranscriptomics-m6a-peak-calling` — `epitranscriptomics/m6a-peak-calling` — usage 41 (exomePeak2, in 1 workflow)
- `bio-epitranscriptomics-m6anet-analysis` — `epitranscriptomics/m6anet-analysis` — usage 22 (m6Anet)

### expression-matrix

- `bio-expression-matrix-metadata-joins` — `expression-matrix/metadata-joins` — usage 50 (pandas)
- `bio-expression-matrix-sparse-handling` — `expression-matrix/sparse-handling` — usage 50 (scipy.sparse)
- `bio-expression-matrix-normalization` — `expression-matrix/normalization` — usage 37 (DESeq2)
- `bio-expression-matrix-gene-id-mapping` — `expression-matrix/gene-id-mapping` — usage 35 (biomaRt)
- `bio-expression-matrix-counts-ingest` — `expression-matrix/counts-ingest` — usage 33 (tximport)

### flow-cytometry

- `bio-flow-cytometry-compensation-transformation` — `flow-cytometry/compensation-transformation` — usage 47 (flowCore, in 1 workflow)
- `bio-flow-cytometry-doublet-detection` — `flow-cytometry/doublet-detection` — usage 47 (flowCore, in 1 workflow)
- `bio-flow-cytometry-fcs-handling` — `flow-cytometry/fcs-handling` — usage 47 (flowCore, in 1 workflow)
- `bio-flow-cytometry-cytometry-qc` — `flow-cytometry/cytometry-qc` — usage 45 (flowAI, in 1 workflow)
- `bio-flow-cytometry-gating-analysis` — `flow-cytometry/gating-analysis` — usage 45 (flowWorkspace, in 1 workflow)
- `bio-flow-cytometry-bead-normalization` — `flow-cytometry/bead-normalization` — usage 43 (CATALYST, in 1 workflow)
- `bio-flow-cytometry-clustering-phenotyping` — `flow-cytometry/clustering-phenotyping` — usage 43 (CATALYST, in 1 workflow)
- `bio-flow-cytometry-differential-analysis` — `flow-cytometry/differential-analysis` — usage 43 (diffcyt, in 1 workflow)

### gene-regulatory-networks

- `bio-gene-regulatory-networks-scenic-regulons` — `gene-regulatory-networks/scenic-regulons` — usage 38 (pySCENIC, in 1 workflow)
- `bio-gene-regulatory-networks-multiomics-grn` — `gene-regulatory-networks/multiomics-grn` — usage 31 (SCENIC+, in 1 workflow)
- `bio-gene-regulatory-networks-coexpression-networks` — `gene-regulatory-networks/coexpression-networks` — usage 31 (WGCNA)
- `bio-gene-regulatory-networks-grn-inference` — `gene-regulatory-networks/grn-inference` — usage 29 (VIPER)
- `bio-gene-regulatory-networks-perturbation-simulation` — `gene-regulatory-networks/perturbation-simulation` — usage 15 (CellOracle, in 1 workflow)
- `bio-gene-regulatory-networks-differential-networks` — `gene-regulatory-networks/differential-networks` — usage 0 (DiffCorr)

### genome-annotation

- `bio-genome-annotation-annotation-qc` — `genome-annotation/annotation-qc` — usage 50 (BUSCO, in 1 workflow)
- `bio-genome-annotation-ncrna-annotation` — `genome-annotation/ncrna-annotation` — usage 50 (Infernal, in 1 workflow)
- `bio-genome-annotation-repeat-annotation` — `genome-annotation/repeat-annotation` — usage 48 (RepeatMasker, in 1 workflow)
- `bio-genome-annotation-functional-annotation` — `genome-annotation/functional-annotation` — usage 47 (eggNOG-mapper, in 1 workflow)
- `bio-genome-annotation-prokaryotic-annotation` — `genome-annotation/prokaryotic-annotation` — usage 47 (Bakta, in 1 workflow)
- `bio-genome-annotation-eukaryotic-gene-prediction` — `genome-annotation/eukaryotic-gene-prediction` — usage 40 (BRAKER3, in 1 workflow)
- `bio-genome-annotation-annotation-transfer` — `genome-annotation/annotation-transfer` — usage 29 (Liftoff)

### genome-assembly

- `bio-genome-assembly-assembly-qc` — `genome-assembly/assembly-qc` — usage 59 (QUAST, in 2 workflows)
- `bio-genome-assembly-short-read-assembly` — `genome-assembly/short-read-assembly` — usage 52 (SPAdes, in 1 workflow)
- `bio-genome-assembly-long-read-assembly` — `genome-assembly/long-read-assembly` — usage 49 (Flye, in 1 workflow)
- `bio-genome-assembly-metagenome-assembly` — `genome-assembly/metagenome-assembly` — usage 49 (metaFlye, in 1 workflow)
- `bio-genome-assembly-assembly-polishing` — `genome-assembly/assembly-polishing` — usage 47 (Pilon, in 1 workflow)
- `bio-genome-assembly-hifi-assembly` — `genome-assembly/hifi-assembly` — usage 47 (hifiasm, in 1 workflow)
- `bio-genome-assembly-genome-profiling` — `genome-assembly/genome-profiling` — usage 43 (GenomeScope2, in 1 workflow)
- `bio-genome-assembly-contamination-detection` — `genome-assembly/contamination-detection` — usage 43 (CheckM2, in 1 workflow)
- `bio-genome-assembly-scaffolding` — `genome-assembly/scaffolding` — usage 42 (YaHS, in 1 workflow)

### genome-engineering

- `bio-genome-engineering-base-editing-design` — `genome-engineering/base-editing-design` — usage 58 (BioPython, in 1 workflow)
- `bio-genome-engineering-hdr-template-design` — `genome-engineering/hdr-template-design` — usage 48 (primer3-py, in 1 workflow)
- `bio-genome-engineering-off-target-prediction` — `genome-engineering/off-target-prediction` — usage 41 (Cas-OFFinder, in 1 workflow)
- `bio-genome-engineering-grna-design` — `genome-engineering/grna-design` — usage 15 (CRISPOR, in 1 workflow)
- `bio-genome-engineering-prime-editing-design` — `genome-engineering/prime-editing-design` — usage 15 (PrimeDesign, in 1 workflow)

### genome-intervals

- `bio-genome-intervals-bed-file-basics` — `genome-intervals/bed-file-basics` — usage 41 (bedtools)
- `bio-genome-intervals-coverage-analysis` — `genome-intervals/coverage-analysis` — usage 41 (bedtools)
- `bio-genome-intervals-interval-arithmetic` — `genome-intervals/interval-arithmetic` — usage 41 (bedtools)
- `bio-genome-intervals-proximity-operations` — `genome-intervals/proximity-operations` — usage 41 (bedtools)
- `bio-genome-intervals-bedgraph-handling` — `genome-intervals/bedgraph-handling` — usage 38 (deeptools)
- `bio-genome-intervals-bigwig-tracks` — `genome-intervals/bigwig-tracks` — usage 38 (pyBigWig)
- `bio-genome-intervals-gtf-gff-handling` — `genome-intervals/gtf-gff-handling` — usage 34 (gffutils)
- `bio-genome-intervals-overlap-significance` — `genome-intervals/overlap-significance` — usage 31 (regioneR)

### hi-c-analysis

- `bio-hi-c-analysis-hic-visualization` — `hi-c-analysis/hic-visualization` — usage 63 (matplotlib, in 1 workflow)
- `bio-hi-c-analysis-compartment-analysis` — `hi-c-analysis/compartment-analysis` — usage 50 (cooltools, in 1 workflow)
- `bio-hi-c-analysis-hic-differential` — `hi-c-analysis/hic-differential` — usage 50 (cooltools, in 1 workflow)
- `bio-hi-c-analysis-loop-calling` — `hi-c-analysis/loop-calling` — usage 50 (cooltools, in 1 workflow)
- `bio-hi-c-analysis-tad-detection` — `hi-c-analysis/tad-detection` — usage 50 (cooltools, in 1 workflow)
- `bio-hi-c-analysis-hic-data-io` — `hi-c-analysis/hic-data-io` — usage 50 (cooler, in 1 workflow)
- `bio-hi-c-analysis-matrix-operations` — `hi-c-analysis/matrix-operations` — usage 50 (cooler, in 1 workflow)
- `bio-hi-c-analysis-contact-pairs` — `hi-c-analysis/contact-pairs` — usage 45 (pairtools, in 1 workflow)
- `bio-hi-c-analysis-hichip-plac-loops` — `hi-c-analysis/hichip-plac-loops` — usage 15 (fithichip, in 1 workflow)

### imaging-mass-cytometry

- `bio-imaging-mass-cytometry-phenotyping` — `imaging-mass-cytometry/phenotyping` — usage 51 (scanpy, in 1 workflow)
- `bio-imaging-mass-cytometry-interactive-annotation` — `imaging-mass-cytometry/interactive-annotation` — usage 49 (napari, in 1 workflow)
- `bio-imaging-mass-cytometry-spatial-analysis` — `imaging-mass-cytometry/spatial-analysis` — usage 45 (squidpy, in 1 workflow)
- `bio-imaging-mass-cytometry-quality-metrics` — `imaging-mass-cytometry/quality-metrics` — usage 43 (CATALYST, in 1 workflow)
- `bio-imaging-mass-cytometry-differential-analysis` — `imaging-mass-cytometry/differential-analysis` — usage 43 (diffcyt, in 1 workflow)
- `bio-imaging-mass-cytometry-data-preprocessing` — `imaging-mass-cytometry/data-preprocessing` — usage 36 (steinbock, in 1 workflow)
- `bio-imaging-mass-cytometry-cell-segmentation` — `imaging-mass-cytometry/cell-segmentation` — usage 15 (deepcell, in 1 workflow)

### immunoinformatics

- `bio-immunoinformatics-mhc-binding-prediction` — `immunoinformatics/mhc-binding-prediction` — usage 44 (mhcflurry, in 1 workflow)
- `bio-immunoinformatics-immunogenicity-scoring` — `immunoinformatics/immunogenicity-scoring` — usage 37 (NeoFox, in 1 workflow)
- `bio-immunoinformatics-neoantigen-prediction` — `immunoinformatics/neoantigen-prediction` — usage 35 (pVACtools, in 1 workflow)
- `bio-immunoinformatics-epitope-prediction` — `immunoinformatics/epitope-prediction` — usage 15 (BepiPred, in 1 workflow)
- `bio-immunoinformatics-mhc-class-ii-prediction` — `immunoinformatics/mhc-class-ii-prediction` — usage 15 (NetMHCIIpan, in 1 workflow)
- `bio-immunoinformatics-tcr-epitope-binding` — `immunoinformatics/tcr-epitope-binding` — usage 0 (tcrdist3)

### liquid-biopsy

- `bio-longitudinal-monitoring` — `liquid-biopsy/longitudinal-monitoring` — usage 65 (pandas, in 1 workflow)
- `bio-ctdna-mutation-detection` — `liquid-biopsy/ctdna-mutation-detection` — usage 47 (VarDict, in 1 workflow)
- `bio-methylation-based-detection` — `liquid-biopsy/methylation-based-detection` — usage 45 (MethylDackel, in 1 workflow)
- `bio-tumor-fraction-estimation` — `liquid-biopsy/tumor-fraction-estimation` — usage 43 (ichorCNA, in 1 workflow)
- `bio-fragment-analysis` — `liquid-biopsy/fragment-analysis` — usage 35 (FinaleToolkit, in 1 workflow)

### long-read-sequencing

- `bio-long-read-sequencing-long-read-qc` — `long-read-sequencing/long-read-qc` — usage 58 (nanoplot, in 2 workflows)
- `bio-long-read-sequencing-long-read-alignment` — `long-read-sequencing/long-read-alignment` — usage 54 (minimap2, in 1 workflow)
- `bio-long-read-sequencing-structural-variants` — `long-read-sequencing/structural-variants` — usage 46 (sniffles, in 1 workflow)
- `bio-long-read-sequencing-medaka-polishing` — `long-read-sequencing/medaka-polishing` — usage 35 (medaka)
- `bio-long-read-sequencing-haplotype-phasing` — `long-read-sequencing/haplotype-phasing` — usage 34 (whatshap)
- `bio-long-read-sequencing-clair3-variants` — `long-read-sequencing/clair3-variants` — usage 31 (Clair3)
- `bio-long-read-sequencing-nanopore-methylation` — `long-read-sequencing/nanopore-methylation` — usage 29 (modkit)
- `bio-long-read-sequencing-isoseq-analysis` — `long-read-sequencing/isoseq-analysis` — usage 17 (SQANTI3)
- `bio-long-read-sequencing-basecalling` — `long-read-sequencing/basecalling` — usage 15 (dorado, in 1 workflow)

### metagenomics

- `bio-metagenomics-metaphlan` — `metagenomics/metaphlan-profiling` — usage 50 (MetaPhlAn, in 1 workflow)
- `bio-metagenomics-visualization` — `metagenomics/metagenome-visualization` — usage 49 (phyloseq, in 1 workflow)
- `bio-metagenomics-contamination-controls` — `metagenomics/contamination-controls` — usage 49 (decontam, in 1 workflow)
- `bio-metagenomics-kraken` — `metagenomics/kraken-classification` — usage 49 (Kraken2, in 1 workflow)
- `bio-metagenomics-functional-profiling` — `metagenomics/functional-profiling` — usage 48 (HUMAnN, in 1 workflow)
- `bio-metagenomics-abundance` — `metagenomics/abundance-estimation` — usage 47 (Bracken, in 1 workflow)
- `bio-metagenomics-amr-detection` — `metagenomics/amr-detection` — usage 33 (AMRFinderPlus)
- `bio-metagenomics-strain-tracking` — `metagenomics/strain-tracking` — usage 31 (inStrain)

### methylation-analysis

- `bio-methylation-differential-cpg` — `methylation-analysis/differential-cpg-testing` — usage 52 (limma, in 1 workflow)
- `bio-methylation-bismark-alignment` — `methylation-analysis/bismark-alignment` — usage 47 (Bismark, in 1 workflow)
- `bio-methylation-calling` — `methylation-analysis/methylation-calling` — usage 47 (Bismark, in 1 workflow)
- `bio-methylation-methylkit` — `methylation-analysis/methylkit-analysis` — usage 46 (methylKit, in 1 workflow)
- `bio-methylation-dmr-detection` — `methylation-analysis/dmr-detection` — usage 43 (dmrseq, in 1 workflow)
- `bio-methylation-array-qc-filtering` — `methylation-analysis/array-qc-filtering` — usage 31 (minfi)
- `bio-methylation-array-preprocessing` — `methylation-analysis/array-preprocessing` — usage 28 (sesame)
- `bio-methylation-cell-type-deconvolution` — `methylation-analysis/cell-type-deconvolution` — usage 28 (EpiDISH)
- `bio-methylation-epigenetic-clocks` — `methylation-analysis/epigenetic-clocks` — usage 26 (methylclock)
- `bio-methylation-ewas-design` — `methylation-analysis/ewas-design` — usage 0 (meffil)

### multi-omics-integration

- `bio-multi-omics-data-harmonization` — `multi-omics-integration/data-harmonization` — usage 49 (MultiAssayExperiment, in 1 workflow)
- `bio-multi-omics-integration-design` — `multi-omics-integration/integration-design` — usage 49 (MultiAssayExperiment, in 1 workflow)
- `bio-multi-omics-similarity-network` — `multi-omics-integration/similarity-network` — usage 47 (SNFtool, in 1 workflow)
- `bio-multi-omics-mixomics-analysis` — `multi-omics-integration/mixomics-analysis` — usage 47 (mixOmics, in 1 workflow)
- `bio-multi-omics-mofa-integration` — `multi-omics-integration/mofa-integration` — usage 42 (MOFA2, in 1 workflow)

### phasing-imputation

- `bio-phasing-imputation-imputation-qc` — `phasing-imputation/imputation-qc` — usage 57 (bcftools, in 1 workflow)
- `bio-phasing-imputation-reference-panels` — `phasing-imputation/reference-panels` — usage 57 (bcftools, in 1 workflow)
- `bio-phasing-imputation-genotype-imputation` — `phasing-imputation/genotype-imputation` — usage 44 (Beagle, in 1 workflow)
- `bio-phasing-imputation-haplotype-phasing` — `phasing-imputation/haplotype-phasing` — usage 40 (SHAPEIT5, in 1 workflow)

### population-genetics

- `bio-population-genetics-plink-basics` — `population-genetics/plink-basics` — usage 48 (plink, in 1 workflow)
- `bio-population-genetics-association-testing` — `population-genetics/association-testing` — usage 47 (plink2, in 1 workflow)
- `bio-population-genetics-linkage-disequilibrium` — `population-genetics/linkage-disequilibrium` — usage 47 (plink2, in 1 workflow)
- `bio-population-genetics-population-structure` — `population-genetics/population-structure` — usage 47 (plink2, in 1 workflow)
- `bio-population-genetics-scikit-allel-analysis` — `population-genetics/scikit-allel-analysis` — usage 40 (scikit-allel)
- `bio-population-genetics-selection-statistics` — `population-genetics/selection-statistics` — usage 40 (scikit-allel)

### primer-design

- `bio-primer-design-primer-basics` — `primer-design/primer-basics` — usage 33 (primer3-py)
- `bio-primer-design-primer-validation` — `primer-design/primer-validation` — usage 33 (primer3-py)
- `bio-primer-design-qpcr-primers` — `primer-design/qpcr-primers` — usage 33 (primer3-py)
- `bio-primer-design-primer-specificity` — `primer-design/primer-specificity` — usage 0 (mfeprimer)

### read-alignment

- `bio-read-alignment-star-alignment` — `read-alignment/star-alignment` — usage 69 (STAR, in 3 workflows)
- `bio-read-alignment-bowtie2-alignment` — `read-alignment/bowtie2-alignment` — usage 65 (bowtie2, in 2 workflows)
- `bio-read-alignment-bwa-alignment` — `read-alignment/bwa-alignment` — usage 56 (bwa-mem2, in 2 workflows)
- `bio-read-alignment-hisat2-alignment` — `read-alignment/hisat2-alignment` — usage 35 (HISAT2)

### read-qc

- `bio-read-qc-fastp-workflow` — `read-qc/fastp-workflow` — usage 87 (fastp, in 9 workflows)
- `bio-read-qc-adapter-trimming` — `read-qc/adapter-trimming` — usage 54 (cutadapt, in 1 workflow)
- `bio-read-qc-quality-reports` — `read-qc/quality-reports` — usage 54 (fastqc, in 1 workflow)
- `bio-read-qc-rnaseq-qc` — `read-qc/rnaseq-qc` — usage 48 (RSeQC, in 1 workflow)
- `bio-read-qc-quality-filtering` — `read-qc/quality-filtering` — usage 35 (trimmomatic)
- `bio-read-qc-umi-processing` — `read-qc/umi-processing` — usage 34 (umi_tools)
- `bio-read-qc-contamination-screening` — `read-qc/contamination-screening` — usage 31 (fastq_screen)

### reporting

- `bio-reporting-figure-export` — `reporting/figure-export` — usage 48 (matplotlib)
- `bio-reporting-rmarkdown-reports` — `reporting/rmarkdown-reports` — usage 41 (rmarkdown)
- `bio-reporting-jupyter-reports` — `reporting/jupyter-reports` — usage 38 (papermill)
- `bio-reporting-automated-qc-reports` — `reporting/automated-qc-reports` — usage 37 (multiqc)
- `bio-reporting-quarto-reports` — `reporting/quarto-reports` — usage 34 (Quarto)
- `bio-reporting-publication-tables` — `reporting/publication-tables` — usage 33 (gtsummary)

### restriction-analysis

- `bio-restriction-enzyme-selection` — `restriction-analysis/enzyme-selection` — usage 43 (Bio.Restriction)
- `bio-restriction-fragment-analysis` — `restriction-analysis/fragment-analysis` — usage 43 (Bio.Restriction)
- `bio-restriction-golden-gate-assembly` — `restriction-analysis/golden-gate-assembly` — usage 43 (Bio.Restriction)
- `bio-restriction-mapping` — `restriction-analysis/restriction-mapping` — usage 43 (Bio.Restriction)
- `bio-restriction-sites` — `restriction-analysis/restriction-sites` — usage 43 (Bio.Restriction)

### ribo-seq

- `bio-ribo-seq-riboseq-preprocessing` — `ribo-seq/riboseq-preprocessing` — usage 54 (STAR, in 1 workflow)
- `bio-ribo-seq-orf-detection` — `ribo-seq/orf-detection` — usage 44 (RiboCode, in 1 workflow)
- `bio-ribo-seq-ribosome-stalling` — `ribo-seq/ribosome-stalling` — usage 44 (Plastid, in 1 workflow)
- `bio-ribo-seq-initiation-site-mapping` — `ribo-seq/initiation-site-mapping` — usage 43 (Ribo-TISH, in 1 workflow)
- `bio-ribo-seq-translation-efficiency` — `ribo-seq/translation-efficiency` — usage 40 (riborex, in 1 workflow)
- `bio-ribo-seq-ribosome-periodicity` — `ribo-seq/ribosome-periodicity` — usage 38 (riboWaltz, in 1 workflow)

### rna-quantification

- `bio-rna-quantification-alignment-free-quant` — `rna-quantification/alignment-free-quant` — usage 61 (salmon, in 2 workflows)
- `bio-rna-quantification-tximport-workflow` — `rna-quantification/tximport-workflow` — usage 57 (tximport, in 2 workflows)
- `bio-rna-quantification-count-matrix-qc` — `rna-quantification/count-matrix-qc` — usage 52 (DESeq2, in 1 workflow)
- `bio-rna-quantification-featurecounts-counting` — `rna-quantification/featurecounts-counting` — usage 35 (featureCounts)

### rna-structure

- `bio-rna-structure-secondary-structure-prediction` — `rna-structure/secondary-structure-prediction` — usage 39 (ViennaRNA)
- `bio-rna-structure-ncrna-search` — `rna-structure/ncrna-search` — usage 35 (Infernal)
- `bio-rna-structure-covariation-analysis` — `rna-structure/covariation-analysis` — usage 32 (R-scape)
- `bio-rna-structure-structure-probing` — `rna-structure/structure-probing` — usage 24 (ShapeMapper2)

### sequence-io

- `bio-compressed-files` — `sequence-io/compressed-files` — usage 43 (Bio.bgzf)
- `bio-fastq-quality` — `sequence-io/fastq-quality` — usage 43 (Bio.SeqIO)
- `bio-filter-sequences` — `sequence-io/filter-sequences` — usage 43 (Bio.SeqIO)
- `bio-format-conversion` — `sequence-io/format-conversion` — usage 43 (Bio.SeqIO)
- `bio-paired-end-fastq` — `sequence-io/paired-end-fastq` — usage 43 (Bio.SeqIO)
- `bio-read-sequences` — `sequence-io/read-sequences` — usage 43 (Bio.SeqIO)
- `bio-sequence-statistics` — `sequence-io/sequence-statistics` — usage 43 (Bio.SeqIO)
- `bio-write-sequences` — `sequence-io/write-sequences` — usage 43 (Bio.SeqIO)

### sequence-manipulation

- `bio-motif-search` — `sequence-manipulation/motif-search` — usage 43 (Bio.motifs)
- `bio-reverse-complement` — `sequence-manipulation/reverse-complement` — usage 43 (Bio.Seq)
- `bio-seq-objects` — `sequence-manipulation/seq-objects` — usage 43 (Bio.Seq)
- `bio-sequence-properties` — `sequence-manipulation/sequence-properties` — usage 43 (Bio.SeqUtils)
- `bio-sequence-slicing` — `sequence-manipulation/sequence-slicing` — usage 43 (Bio.Seq)
- `bio-transcription-translation` — `sequence-manipulation/transcription-translation` — usage 43 (Bio.Seq)

### small-rna-seq

- `bio-small-rna-seq-smrna-preprocessing` — `small-rna-seq/smrna-preprocessing` — usage 54 (cutadapt, in 1 workflow)
- `bio-small-rna-seq-differential-mirna` — `small-rna-seq/differential-mirna` — usage 52 (DESeq2, in 1 workflow)
- `bio-small-rna-seq-mirdeep2-analysis` — `small-rna-seq/mirdeep2-analysis` — usage 44 (miRDeep2, in 1 workflow)
- `bio-small-rna-seq-mirge3-analysis` — `small-rna-seq/mirge3-analysis` — usage 43 (miRge3, in 1 workflow)
- `bio-small-rna-seq-target-prediction` — `small-rna-seq/target-prediction` — usage 42 (miRanda, in 1 workflow)
- `bio-small-rna-seq-trf-pirna-profiling` — `small-rna-seq/trf-pirna-profiling` — usage 40 (MINTmap, in 1 workflow)

### spatial-transcriptomics

- `bio-spatial-transcriptomics-image-analysis` — `spatial-transcriptomics/image-analysis` — usage 45 (squidpy, in 1 workflow)
- `bio-spatial-transcriptomics-spatial-domains` — `spatial-transcriptomics/spatial-domains` — usage 45 (squidpy, in 1 workflow)
- `bio-spatial-transcriptomics-spatial-neighbors` — `spatial-transcriptomics/spatial-neighbors` — usage 45 (squidpy, in 1 workflow)
- `bio-spatial-transcriptomics-spatial-preprocessing` — `spatial-transcriptomics/spatial-preprocessing` — usage 45 (squidpy, in 1 workflow)
- `bio-spatial-transcriptomics-spatial-statistics` — `spatial-transcriptomics/spatial-statistics` — usage 45 (squidpy, in 1 workflow)
- `bio-spatial-transcriptomics-spatial-visualization` — `spatial-transcriptomics/spatial-visualization` — usage 45 (squidpy, in 1 workflow)
- `bio-spatial-transcriptomics-spatial-data-io` — `spatial-transcriptomics/spatial-data-io` — usage 45 (spatialdata, in 1 workflow)
- `bio-spatial-transcriptomics-spatial-communication` — `spatial-transcriptomics/spatial-communication` — usage 30 (squidpy)
- `bio-spatial-transcriptomics-spatial-multiomics` — `spatial-transcriptomics/spatial-multiomics` — usage 29 (muon)
- `bio-spatial-transcriptomics-high-resolution-binning` — `spatial-transcriptomics/high-resolution-binning` — usage 19 (bin2cell)
- `bio-spatial-transcriptomics-spatial-deconvolution` — `spatial-transcriptomics/spatial-deconvolution` — usage 15 (cell2location, in 1 workflow)
- `bio-spatial-transcriptomics-spatial-proteomics` — `spatial-transcriptomics/spatial-proteomics` — usage 0 (scimap)

### structural-biology

- `bio-structural-biology-alphafold-predictions` — `structural-biology/alphafold-predictions` — usage 50 (requests)
- `bio-structural-biology-geometric-analysis` — `structural-biology/geometric-analysis` — usage 43 (Bio.PDB)
- `bio-structural-biology-interface-analysis` — `structural-biology/interface-analysis` — usage 43 (Bio.PDB)
- `bio-structural-biology-structure-io` — `structural-biology/structure-io` — usage 43 (Bio.PDB)
- `bio-structural-biology-structure-modification` — `structural-biology/structure-modification` — usage 43 (Bio.PDB)
- `bio-structural-biology-structure-navigation` — `structural-biology/structure-navigation` — usage 43 (Bio.PDB)
- `bio-structural-biology-structure-validation` — `structural-biology/structure-validation` — usage 43 (Bio.PDB)
- `bio-structural-biology-structure-preparation` — `structural-biology/structure-preparation` — usage 39 (PDBFixer)
- `bio-structural-biology-binding-site-detection` — `structural-biology/binding-site-detection` — usage 32 (fpocket)
- `bio-structural-biology-modern-structure-prediction` — `structural-biology/modern-structure-prediction` — usage 0 (ESMFold)

### systems-biology

- `bio-systems-biology-context-specific-models` — `systems-biology/context-specific-models` — usage 46 (cobrapy, in 1 workflow)
- `bio-systems-biology-flux-balance-analysis` — `systems-biology/flux-balance-analysis` — usage 46 (cobrapy, in 1 workflow)
- `bio-systems-biology-gene-essentiality` — `systems-biology/gene-essentiality` — usage 46 (cobrapy, in 1 workflow)
- `bio-systems-biology-metabolic-reconstruction` — `systems-biology/metabolic-reconstruction` — usage 38 (CarveMe, in 1 workflow)
- `bio-systems-biology-model-curation` — `systems-biology/model-curation` — usage 36 (memote, in 1 workflow)
- `bio-systems-biology-community-metabolic-modeling` — `systems-biology/community-metabolic-modeling` — usage 29 (micom)
- `bio-systems-biology-strain-design` — `systems-biology/strain-design` — usage 25 (straindesign)

### tcr-bcr-analysis

- `bio-tcr-bcr-analysis-immcantation-analysis` — `tcr-bcr-analysis/immcantation-analysis` — usage 48 (alakazam, in 1 workflow)
- `bio-tcr-bcr-analysis-scirpy-analysis` — `tcr-bcr-analysis/scirpy-analysis` — usage 46 (scirpy, in 1 workflow)
- `bio-tcr-bcr-analysis-mixcr-analysis` — `tcr-bcr-analysis/mixcr-analysis` — usage 15 (MiXCR, in 1 workflow)
- `bio-tcr-bcr-analysis-repertoire-visualization` — `tcr-bcr-analysis/repertoire-visualization` — usage 15 (VDJtools, in 1 workflow)
- `bio-tcr-bcr-analysis-specificity-annotation` — `tcr-bcr-analysis/specificity-annotation` — usage 15 (tcrdist3, in 1 workflow)
- `bio-tcr-bcr-analysis-vdjtools-analysis` — `tcr-bcr-analysis/vdjtools-analysis` — usage 15 (VDJtools, in 1 workflow)

### temporal-genomics

- `bio-temporal-genomics-trajectory-modeling` — `temporal-genomics/trajectory-modeling` — usage 57 (mgcv, in 1 workflow)
- `bio-temporal-genomics-periodicity-detection` — `temporal-genomics/periodicity-detection` — usage 50 (scipy)
- `bio-temporal-genomics-temporal-grn` — `temporal-genomics/temporal-grn` — usage 47 (statsmodels)
- `bio-temporal-genomics-temporal-clustering` — `temporal-genomics/temporal-clustering` — usage 43 (Mfuzz, in 1 workflow)
- `bio-temporal-genomics-circadian-rhythms` — `temporal-genomics/circadian-rhythms` — usage 15 (CosinorPy, in 1 workflow)
- `bio-temporal-genomics-differential-rhythmicity` — `temporal-genomics/differential-rhythmicity` — usage 0 (limorhyde)

### variant-calling

- `bio-variant-calling` — `variant-calling/variant-calling` — usage 57 (bcftools, in 1 workflow)
- `bio-variant-calling-clinical-interpretation` — `variant-calling/clinical-interpretation` — usage 57 (bcftools, in 1 workflow)
- `bio-gatk-variant-calling` — `variant-calling/gatk-variant-calling` — usage 53 (gatk, in 1 workflow)
- `bio-variant-calling-joint-calling` — `variant-calling/joint-calling` — usage 53 (GATK, in 1 workflow)
- `bio-variant-calling-structural-variant-calling` — `variant-calling/structural-variant-calling` — usage 47 (manta, in 1 workflow)
- `bio-consensus-sequences` — `variant-calling/consensus-sequences` — usage 42 (bcftools)
- `bio-variant-calling-deepvariant` — `variant-calling/deepvariant` — usage 31 (DeepVariant)

### workflow-management

- `bio-workflow-management-cwl-workflows` — `workflow-management/cwl-workflows` — usage 40 (cwltool)
- `bio-workflow-management-snakemake-workflows` — `workflow-management/snakemake-workflows` — usage 39 (Snakemake)
- `bio-workflow-management-nextflow-pipelines` — `workflow-management/nextflow-pipelines` — usage 36 (Nextflow)
- `bio-workflow-management-wdl-workflows` — `workflow-management/wdl-workflows` — usage 34 (cromwell)
- `bio-workflow-management-nf-core-pipelines` — `workflow-management/nf-core-pipelines` — usage 32 (nf-core)

### workflows

- `bio-workflows-biomarker-pipeline` — `workflows/biomarker-pipeline` — usage 48 (sklearn)
- `bio-workflows-clinical-trial-pipeline` — `workflows/clinical-trial-pipeline` — usage 47 (statsmodels)
- `bio-workflows-fastq-to-variants` — `workflows/fastq-to-variants` — usage 42 (bcftools)
- `bio-workflows-riboseq-pipeline` — `workflows/riboseq-pipeline` — usage 39 (STAR)
- `bio-workflows-somatic-variant-pipeline` — `workflows/somatic-variant-pipeline` — usage 38 (GATK Mutect2)
- `bio-workflows-rnaseq-to-de` — `workflows/rnaseq-to-de` — usage 37 (DESeq2)
- `bio-workflows-microbiome-pipeline` — `workflows/microbiome-pipeline` — usage 37 (DADA2)
- `bio-workflows-multiome-pipeline` — `workflows/multiome-pipeline` — usage 36 (Seurat)
- `bio-workflows-outbreak-pipeline` — `workflows/outbreak-pipeline` — usage 35 (mlst)
- `bio-workflows-hic-pipeline` — `workflows/hic-pipeline` — usage 34 (cooler)
- `bio-workflows-cnv-pipeline` — `workflows/cnv-pipeline` — usage 34 (CNVkit)
- `bio-workflows-genome-assembly-pipeline` — `workflows/genome-assembly-pipeline` — usage 34 (Flye)
- `bio-workflows-metagenomics-pipeline` — `workflows/metagenomics-pipeline` — usage 34 (Kraken2)
- `bio-workflows-expression-to-pathways` — `workflows/expression-to-pathways` — usage 33 (clusterProfiler)
- `bio-workflows-genome-annotation-pipeline` — `workflows/genome-annotation-pipeline` — usage 32 (Bakta)
- `bio-workflows-methylation-pipeline` — `workflows/methylation-pipeline` — usage 32 (Bismark)
- `bio-workflows-gwas-pipeline` — `workflows/gwas-pipeline` — usage 32 (PLINK2)
- `bio-workflows-longread-sv-pipeline` — `workflows/longread-sv-pipeline` — usage 31 (Sniffles)
- `bio-workflows-splicing-pipeline` — `workflows/splicing-pipeline` — usage 31 (rMATS-turbo)
- `bio-workflows-metabolic-modeling-pipeline` — `workflows/metabolic-modeling-pipeline` — usage 30 (cobrapy)
- `bio-workflows-spatial-pipeline` — `workflows/spatial-pipeline` — usage 30 (Squidpy)
- `bio-workflows-atacseq-pipeline` — `workflows/atacseq-pipeline` — usage 29 (MACS3)
- `bio-workflows-chipseq-pipeline` — `workflows/chipseq-pipeline` — usage 29 (MACS3)
- `bio-workflows-liquid-biopsy-pipeline` — `workflows/liquid-biopsy-pipeline` — usage 28 (ichorCNA)
- `bio-workflows-timecourse-pipeline` — `workflows/timecourse-pipeline` — usage 28 (Mfuzz)
- `bio-workflows-smrna-pipeline` — `workflows/smrna-pipeline` — usage 28 (miRge3.0)
- `bio-workflows-cytometry-pipeline` — `workflows/cytometry-pipeline` — usage 28 (CATALYST)
- `bio-workflows-edna-pipeline` — `workflows/edna-pipeline` — usage 28 (obitools3)
- `bio-workflows-multi-omics-pipeline` — `workflows/multi-omics-pipeline` — usage 27 (MOFA2)
- `bio-workflows-merip-pipeline` — `workflows/merip-pipeline` — usage 26 (exomePeak2)
- `bio-workflows-grn-pipeline` — `workflows/grn-pipeline` — usage 23 (pySCENIC)
- `bio-workflows-clip-pipeline` — `workflows/clip-pipeline` — usage 21 (CLIPper)
- `bio-workflows-imc-pipeline` — `workflows/imc-pipeline` — usage 21 (steinbock)
- `bio-workflows-neoantigen-pipeline` — `workflows/neoantigen-pipeline` — usage 20 (pVACtools)
- `bio-workflows-causal-genomics-pipeline` — `workflows/causal-genomics-pipeline` — usage 0 (TwoSampleMR)
- `bio-workflows-crispr-editing-pipeline` — `workflows/crispr-editing-pipeline` — usage 0 (CRISPOR)
- `bio-workflows-tcr-pipeline` — `workflows/tcr-pipeline` — usage 0 (MiXCR)
