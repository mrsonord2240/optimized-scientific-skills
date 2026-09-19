# Remaining Skills

Not yet refined. Scope is deliberately limited to the rest of
[GPTomics/bioSkills](https://github.com/GPTomics/bioSkills) at commit
`d91ed3d563019e649dc854c56ccd62551359488a`; other source corpora are out of scope for now.

**456 remaining** across 59 folders. 105 are already refined and live in `skills/`.

The source tree holds 562 Skills: 105 refined, 0 audited and excluded, 1 out of scope, 456 remaining.

| folder | remaining | refined |
| --- | ---: | ---: |
| workflows | 37 | 4 |
| data-visualization | 20 | 0 |
| comparative-genomics | 13 | 0 |
| atac-seq | 12 | 0 |
| chemoinformatics | 12 | 8 |
| chip-seq | 12 | 0 |
| clinical-biostatistics | 12 | 0 |
| clip-seq | 12 | 0 |
| spatial-transcriptomics | 12 | 0 |
| copy-number | 11 | 0 |
| alignment-files | 10 | 0 |
| methylation-analysis | 10 | 0 |
| single-cell | 10 | 7 |
| structural-biology | 10 | 0 |
| alternative-splicing | 9 | 0 |
| genome-assembly | 9 | 0 |
| hi-c-analysis | 9 | 0 |
| long-read-sequencing | 9 | 0 |
| sequence-io | 9 | 0 |
| clinical-databases | 8 | 4 |
| flow-cytometry | 8 | 0 |
| genome-intervals | 8 | 0 |
| metagenomics | 8 | 0 |
| genome-annotation | 7 | 0 |
| imaging-mass-cytometry | 7 | 0 |
| liquid-biopsy | 7 | 0 |
| read-qc | 7 | 0 |
| sequence-manipulation | 7 | 0 |
| systems-biology | 7 | 0 |
| variant-calling | 7 | 6 |
| database-access | 6 | 9 |
| ecological-genomics | 6 | 0 |
| gene-regulatory-networks | 6 | 0 |
| immunoinformatics | 6 | 0 |
| population-genetics | 6 | 1 |
| reporting | 6 | 0 |
| ribo-seq | 6 | 0 |
| small-rna-seq | 6 | 0 |
| tcr-bcr-analysis | 6 | 0 |
| temporal-genomics | 6 | 0 |
| alignment | 5 | 2 |
| differential-expression | 5 | 1 |
| epidemiological-genomics | 5 | 0 |
| epitranscriptomics | 5 | 0 |
| expression-matrix | 5 | 0 |
| genome-engineering | 5 | 0 |
| multi-omics-integration | 5 | 0 |
| restriction-analysis | 5 | 0 |
| workflow-management | 5 | 0 |
| machine-learning | 4 | 2 |
| phasing-imputation | 4 | 0 |
| primer-design | 4 | 0 |
| read-alignment | 4 | 0 |
| rna-quantification | 4 | 0 |
| rna-structure | 4 | 0 |
| crispr-screens | 3 | 12 |
| microbiome | 3 | 3 |
| causal-genomics | 1 | 10 |
| proteomics | 1 | 8 |

## Promoted, fix pass still needed

These Skills are in `skills/` because their audit found them deployable with no open P0.
They have not yet been through a fix pass, however high they scored. Their open findings
are in the audit record. `fix_pass` in PROVENANCE.json carries the same flag.

| skill | score | grade |
| --- | ---: | --- |
| `bio-admet-prediction` | 86 | Limited Release |
| `bio-causal-genomics-genetic-correlation` | 89 | Production Ready |
| `bio-differential-expression-deseq2-basics` | 92 | Production Ready |
| `bio-machine-learning-model-validation` | 93 | Limited Release |
| `bio-machine-learning-prediction-explanation` | 89 | Limited Release |
| `bio-molecular-descriptors` | 86 | Limited Release |
| `bio-molecular-io` | 91 | Limited Release |
| `bio-molecular-standardization` | 90 | Limited Release |
| `bio-qsar-modeling` | 88 | Limited Release |
| `bio-scaffold-analysis` | 88 | Limited Release |
| `bio-similarity-searching` | 90 | Limited Release |
| `bio-single-cell-batch-integration` | 90 | Limited Release |
| `bio-single-cell-cell-annotation` | 87 | Limited Release |
| `bio-single-cell-clustering` | 89 | Limited Release |
| `bio-single-cell-differential-abundance` | 88 | Production Ready |
| `bio-single-cell-doublet-detection` | 85 | Limited Release |
| `bio-single-cell-markers-annotation` | 88 | Production Ready |
| `bio-single-cell-preprocessing` | 85 | Limited Release |
| `bio-substructure-search` | 88 | Limited Release |
| `bio-workflows-scrnaseq-pipeline` | 85 | Limited Release |

## Promoted, re-audit still needed

Changed in staging after their latest audit, so the score below describes earlier bytes.
`reaudit` in PROVENANCE.json carries the same flag.

| skill | score at last audit | audited on |
| --- | ---: | --- |
| `bio-alignment-multiple` | 85 | 2026-09-11 |
| `bio-alignment-trimming` | 83 | 2026-09-15 |
| `bio-causal-genomics-genetic-correlation` | 89 | 2026-09-17 |
| `bio-clinical-databases-clinvar-lookup` | 84 | 2026-09-15 |
| `bio-clinical-databases-dbsnp-queries` | 86 | 2026-09-15 |
| `bio-clinical-databases-gnomad-frequencies` | 86 | 2026-09-15 |
| `bio-clinical-databases-myvariant-queries` | 84 | 2026-09-15 |
| `bio-phylo-bayesian-inference` | 87 | 2026-09-15 |
| `bio-phylo-distance-calculations` | 87 | 2026-09-15 |
| `bio-phylo-divergence-dating` | 87 | 2026-09-15 |
| `bio-phylo-modern-tree-inference` | 88 | 2026-09-15 |
| `bio-phylo-species-trees` | 86 | 2026-09-15 |
| `bio-phylo-tree-io` | 88 | 2026-09-15 |
| `bio-phylo-tree-manipulation` | 86 | 2026-09-15 |
| `bio-population-genetics-rare-variant-association` | 89 | 2026-09-15 |
| `bio-variant-annotation` | 86 | 2026-09-15 |
| `bio-variant-calling-filtering-best-practices` | 88 | 2026-09-15 |
| `bio-variant-normalization` | 89 | 2026-09-15 |
| `bio-vcf-basics` | 90 | 2026-09-11 |
| `bio-vcf-manipulation` | 90 | 2026-09-15 |
| `bio-vcf-statistics` | 88 | 2026-09-15 |

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

### alignment

- `bio-alignment-io` — `alignment/alignment-io`
- `bio-alignment-msa-parsing` — `alignment/msa-parsing`
- `bio-alignment-msa-statistics` — `alignment/msa-statistics`
- `bio-alignment-pairwise` — `alignment/pairwise-alignment`
- `bio-alignment-structural` — `alignment/structural-alignment`

### alignment-files

- `bio-alignment-amplicon-clipping` — `alignment-files/alignment-amplicon-clipping`
- `bio-alignment-filtering` — `alignment-files/alignment-filtering`
- `bio-alignment-indexing` — `alignment-files/alignment-indexing`
- `bio-alignment-sorting` — `alignment-files/alignment-sorting`
- `bio-alignment-validation` — `alignment-files/alignment-validation`
- `bio-bam-statistics` — `alignment-files/bam-statistics`
- `bio-duplicate-handling` — `alignment-files/duplicate-handling`
- `bio-pileup-generation` — `alignment-files/pileup-generation`
- `bio-reference-operations` — `alignment-files/reference-operations`
- `bio-sam-bam-basics` — `alignment-files/sam-bam-basics`

### alternative-splicing

- `bio-differential-splicing` — `alternative-splicing/differential-splicing`
- `bio-isoform-switching` — `alternative-splicing/isoform-switching`
- `bio-long-read-splicing` — `alternative-splicing/long-read-splicing`
- `bio-outlier-splicing-detection` — `alternative-splicing/outlier-splicing-detection`
- `bio-sashimi-plots` — `alternative-splicing/sashimi-plots`
- `bio-single-cell-splicing` — `alternative-splicing/single-cell-splicing`
- `bio-splice-variant-prediction` — `alternative-splicing/splice-variant-prediction`
- `bio-splicing-qc` — `alternative-splicing/splicing-qc`
- `bio-splicing-quantification` — `alternative-splicing/splicing-quantification`

### atac-seq

- `bio-atac-seq-allele-specific-accessibility` — `atac-seq/allele-specific-accessibility`
- `bio-atac-seq-atac-peak-calling` — `atac-seq/atac-peak-calling`
- `bio-atac-seq-atac-qc` — `atac-seq/atac-qc`
- `bio-atac-seq-co-accessibility` — `atac-seq/co-accessibility`
- `bio-atac-seq-consensus-peakset` — `atac-seq/consensus-peakset`
- `bio-atac-seq-deep-learning-atac` — `atac-seq/deep-learning-atac`
- `bio-atac-seq-differential-accessibility` — `atac-seq/differential-accessibility`
- `bio-atac-seq-enhancer-gene-linking` — `atac-seq/enhancer-gene-linking`
- `bio-atac-seq-footprinting` — `atac-seq/footprinting`
- `bio-atac-seq-motif-deviation` — `atac-seq/motif-deviation`
- `bio-atac-seq-nucleosome-positioning` — `atac-seq/nucleosome-positioning`
- `bio-atac-seq-single-cell-atac` — `atac-seq/single-cell-atac`

### causal-genomics

- `bio-causal-genomics-transcriptome-wide-association` — `causal-genomics/transcriptome-wide-association`

### chemoinformatics

- `bio-conformer-generation` — `chemoinformatics/conformer-generation`
- `bio-covalent-design` — `chemoinformatics/covalent-design`
- `bio-free-energy-calculations` — `chemoinformatics/free-energy-calculations`
- `bio-generative-design` — `chemoinformatics/generative-design`
- `bio-ml-docking-rescoring` — `chemoinformatics/ml-docking-rescoring`
- `bio-pharmacophore-modeling` — `chemoinformatics/pharmacophore-modeling`
- `bio-pose-validation` — `chemoinformatics/pose-validation`
- `bio-protac-degraders` — `chemoinformatics/protac-degraders`
- `bio-reaction-enumeration` — `chemoinformatics/reaction-enumeration`
- `bio-retrosynthesis` — `chemoinformatics/retrosynthesis`
- `bio-shape-similarity` — `chemoinformatics/shape-similarity`
- `bio-virtual-screening` — `chemoinformatics/virtual-screening`

### chip-seq

- `bio-chipseq-allele-specific-binding` — `chip-seq/allele-specific-binding`
- `bio-chipseq-chip-deep-learning` — `chip-seq/chip-deep-learning`
- `bio-chipseq-chromatin-state-segmentation` — `chip-seq/chromatin-state-segmentation`
- `bio-chipseq-cut-and-run-tag` — `chip-seq/cut-and-run-tag`
- `bio-chipseq-differential-binding` — `chip-seq/differential-binding`
- `bio-chipseq-motif-analysis` — `chip-seq/motif-analysis`
- `bio-chipseq-peak-annotation` — `chip-seq/peak-annotation`
- `bio-chipseq-peak-calling` — `chip-seq/peak-calling`
- `bio-chipseq-qc` — `chip-seq/chipseq-qc`
- `bio-chipseq-spike-in-normalization` — `chip-seq/spike-in-normalization`
- `bio-chipseq-super-enhancers` — `chip-seq/super-enhancers`
- `bio-chipseq-visualization` — `chip-seq/chipseq-visualization`

### clinical-biostatistics

- `bio-clinical-biostatistics-adaptive-designs` — `clinical-biostatistics/adaptive-designs`
- `bio-clinical-biostatistics-bayesian-trials` — `clinical-biostatistics/bayesian-trials`
- `bio-clinical-biostatistics-categorical-tests` — `clinical-biostatistics/categorical-tests`
- `bio-clinical-biostatistics-cdisc-data` — `clinical-biostatistics/cdisc-data-handling`
- `bio-clinical-biostatistics-effect-measures` — `clinical-biostatistics/effect-measures`
- `bio-clinical-biostatistics-logistic-regression` — `clinical-biostatistics/logistic-regression`
- `bio-clinical-biostatistics-missing-data` — `clinical-biostatistics/missing-data-sensitivity`
- `bio-clinical-biostatistics-multiplicity-graphical` — `clinical-biostatistics/multiplicity-graphical`
- `bio-clinical-biostatistics-power-sample-size` — `clinical-biostatistics/power-and-sample-size`
- `bio-clinical-biostatistics-subgroup-analysis` — `clinical-biostatistics/subgroup-analysis`
- `bio-clinical-biostatistics-survival-analysis` — `clinical-biostatistics/survival-analysis`
- `bio-clinical-biostatistics-trial-reporting` — `clinical-biostatistics/trial-reporting`

### clinical-databases

- `bio-clinical-databases-acmg-classification` — `clinical-databases/acmg-classification`
- `bio-clinical-databases-hla-typing` — `clinical-databases/hla-typing`
- `bio-clinical-databases-msi-detection` — `clinical-databases/msi-detection`
- `bio-clinical-databases-pharmacogenomics` — `clinical-databases/pharmacogenomics`
- `bio-clinical-databases-polygenic-risk` — `clinical-databases/polygenic-risk`
- `bio-clinical-databases-somatic-signatures` — `clinical-databases/somatic-signatures`
- `bio-clinical-databases-tumor-mutational-burden` — `clinical-databases/tumor-mutational-burden`
- `bio-clinical-databases-variant-prioritization` — `clinical-databases/variant-prioritization`

### clip-seq

- `bio-clip-seq-ago-clip-mirna-targets` — `clip-seq/ago-clip-mirna-targets`
- `bio-clip-seq-binding-site-annotation` — `clip-seq/binding-site-annotation`
- `bio-clip-seq-clip-alignment` — `clip-seq/clip-alignment`
- `bio-clip-seq-clip-deep-learning` — `clip-seq/clip-deep-learning`
- `bio-clip-seq-clip-motif-analysis` — `clip-seq/clip-motif-analysis`
- `bio-clip-seq-clip-peak-calling` — `clip-seq/clip-peak-calling`
- `bio-clip-seq-clip-preprocessing` — `clip-seq/clip-preprocessing`
- `bio-clip-seq-clip-qc` — `clip-seq/clip-qc`
- `bio-clip-seq-crosslink-site-detection` — `clip-seq/crosslink-site-detection`
- `bio-clip-seq-differential-clip` — `clip-seq/differential-clip`
- `bio-clip-seq-m6a-clip` — `clip-seq/m6a-clip`
- `bio-clip-seq-stamp-antibody-free` — `clip-seq/stamp-antibody-free`

### comparative-genomics

- `bio-comparative-genomics-ancestral-reconstruction` — `comparative-genomics/ancestral-reconstruction`
- `bio-comparative-genomics-comparative-annotation-projection` — `comparative-genomics/comparative-annotation-projection`
- `bio-comparative-genomics-gene-family-evolution` — `comparative-genomics/gene-family-evolution`
- `bio-comparative-genomics-gene-tree-species-tree-reconciliation` — `comparative-genomics/gene-tree-species-tree-reconciliation`
- `bio-comparative-genomics-genome-distance-and-species-delineation` — `comparative-genomics/genome-distance-and-species-delineation`
- `bio-comparative-genomics-hgt-detection` — `comparative-genomics/hgt-detection`
- `bio-comparative-genomics-introgression-detection` — `comparative-genomics/introgression-detection`
- `bio-comparative-genomics-ortholog-inference` — `comparative-genomics/ortholog-inference`
- `bio-comparative-genomics-pangenome-analysis` — `comparative-genomics/pangenome-analysis`
- `bio-comparative-genomics-positive-selection` — `comparative-genomics/positive-selection`
- `bio-comparative-genomics-synteny-analysis` — `comparative-genomics/synteny-analysis`
- `bio-comparative-genomics-whole-genome-alignment` — `comparative-genomics/whole-genome-alignment`
- `bio-comparative-genomics-whole-genome-duplication` — `comparative-genomics/whole-genome-duplication`

### copy-number

- `bio-copy-number-allele-specific-copy-number` — `copy-number/allele-specific-copy-number`
- `bio-copy-number-cnv-annotation` — `copy-number/cnv-annotation`
- `bio-copy-number-cnv-visualization` — `copy-number/cnv-visualization`
- `bio-copy-number-cnvkit-analysis` — `copy-number/cnvkit-analysis`
- `bio-copy-number-copy-ratio-segmentation` — `copy-number/copy-ratio-segmentation`
- `bio-copy-number-focal-amplification-ecdna` — `copy-number/focal-amplification-ecdna`
- `bio-copy-number-gatk-cnv` — `copy-number/gatk-cnv`
- `bio-copy-number-germline-cnv-interpretation` — `copy-number/germline-cnv-interpretation`
- `bio-copy-number-hrd-scoring` — `copy-number/hrd-scoring`
- `bio-copy-number-recurrent-cnv` — `copy-number/recurrent-cnv`
- `bio-copy-number-subclonal-copy-number` — `copy-number/subclonal-copy-number`

### crispr-screens

- `bio-crispr-screens-combinatorial-screens` — `crispr-screens/combinatorial-screens`
- `bio-crispr-screens-in-vivo-screens` — `crispr-screens/in-vivo-screens`
- `bio-crispr-screens-perturb-seq-analysis` — `crispr-screens/perturb-seq-analysis`

### data-visualization

- `bio-data-visualization-circos-plots` — `data-visualization/circos-plots`
- `bio-data-visualization-color-palettes` — `data-visualization/color-palettes`
- `bio-data-visualization-dimensionality-reduction-plots` — `data-visualization/dimensionality-reduction-plots`
- `bio-data-visualization-distribution-plots` — `data-visualization/distribution-plots`
- `bio-data-visualization-flow-and-transition-plots` — `data-visualization/flow-and-transition-plots`
- `bio-data-visualization-forest-funnel-plots` — `data-visualization/forest-funnel-plots`
- `bio-data-visualization-genome-tracks` — `data-visualization/genome-tracks`
- `bio-data-visualization-ggplot2-fundamentals` — `data-visualization/ggplot2-fundamentals`
- `bio-data-visualization-heatmaps-clustering` — `data-visualization/heatmaps-clustering`
- `bio-data-visualization-interactive-visualization` — `data-visualization/interactive-visualization`
- `bio-data-visualization-lollipop-protein-maps` — `data-visualization/lollipop-protein-maps`
- `bio-data-visualization-manhattan-qq-locuszoom` — `data-visualization/manhattan-qq-locuszoom`
- `bio-data-visualization-matplotlib-fundamentals` — `data-visualization/matplotlib-fundamentals`
- `bio-data-visualization-multipanel-figures` — `data-visualization/multipanel-figures`
- `bio-data-visualization-network-visualization` — `data-visualization/network-visualization`
- `bio-data-visualization-oncoprint-mutation-matrices` — `data-visualization/oncoprint-mutation-matrices`
- `bio-data-visualization-sequence-logos` — `data-visualization/sequence-logos`
- `bio-data-visualization-statistical-annotation` — `data-visualization/statistical-annotation`
- `bio-data-visualization-upset-plots` — `data-visualization/upset-plots`
- `bio-data-visualization-volcano-and-ma-plots` — `data-visualization/volcano-and-ma-plots`

### database-access

- `bio-blast-searches` — `database-access/blast-searches`
- `bio-ensembl-rest` — `database-access/ensembl-rest`
- `bio-interaction-databases` — `database-access/interaction-databases`
- `bio-ncbi-datasets-cli` — `database-access/ncbi-datasets-cli`
- `bio-ortholog-inference` — `database-access/ortholog-inference`
- `bio-uniprot-access` — `database-access/uniprot-access`

### differential-expression

- `bio-differential-expression-batch-correction` — `differential-expression/batch-correction`
- `bio-differential-expression-de-results` — `differential-expression/de-results`
- `bio-differential-expression-de-visualization` — `differential-expression/de-visualization`
- `bio-differential-expression-edger-basics` — `differential-expression/edger-basics`
- `bio-differential-expression-timeseries-de` — `differential-expression/timeseries-de`

### ecological-genomics

- `bio-ecological-genomics-biodiversity-metrics` — `ecological-genomics/biodiversity-metrics`
- `bio-ecological-genomics-community-ecology` — `ecological-genomics/community-ecology`
- `bio-ecological-genomics-conservation-genetics` — `ecological-genomics/conservation-genetics`
- `bio-ecological-genomics-edna-metabarcoding` — `ecological-genomics/edna-metabarcoding`
- `bio-ecological-genomics-landscape-genomics` — `ecological-genomics/landscape-genomics`
- `bio-ecological-genomics-species-delimitation` — `ecological-genomics/species-delimitation`

### epidemiological-genomics

- `bio-epidemiological-genomics-amr-surveillance` — `epidemiological-genomics/amr-surveillance`
- `bio-epidemiological-genomics-pathogen-typing` — `epidemiological-genomics/pathogen-typing`
- `bio-epidemiological-genomics-phylodynamics` — `epidemiological-genomics/phylodynamics`
- `bio-epidemiological-genomics-transmission-inference` — `epidemiological-genomics/transmission-inference`
- `bio-epidemiological-genomics-variant-surveillance` — `epidemiological-genomics/variant-surveillance`

### epitranscriptomics

- `bio-epitranscriptomics-m6a-differential` — `epitranscriptomics/m6a-differential`
- `bio-epitranscriptomics-m6a-peak-calling` — `epitranscriptomics/m6a-peak-calling`
- `bio-epitranscriptomics-m6anet-analysis` — `epitranscriptomics/m6anet-analysis`
- `bio-epitranscriptomics-merip-preprocessing` — `epitranscriptomics/merip-preprocessing`
- `bio-epitranscriptomics-modification-visualization` — `epitranscriptomics/modification-visualization`

### expression-matrix

- `bio-expression-matrix-counts-ingest` — `expression-matrix/counts-ingest`
- `bio-expression-matrix-gene-id-mapping` — `expression-matrix/gene-id-mapping`
- `bio-expression-matrix-metadata-joins` — `expression-matrix/metadata-joins`
- `bio-expression-matrix-normalization` — `expression-matrix/normalization`
- `bio-expression-matrix-sparse-handling` — `expression-matrix/sparse-handling`

### flow-cytometry

- `bio-flow-cytometry-bead-normalization` — `flow-cytometry/bead-normalization`
- `bio-flow-cytometry-clustering-phenotyping` — `flow-cytometry/clustering-phenotyping`
- `bio-flow-cytometry-compensation-transformation` — `flow-cytometry/compensation-transformation`
- `bio-flow-cytometry-cytometry-qc` — `flow-cytometry/cytometry-qc`
- `bio-flow-cytometry-differential-analysis` — `flow-cytometry/differential-analysis`
- `bio-flow-cytometry-doublet-detection` — `flow-cytometry/doublet-detection`
- `bio-flow-cytometry-fcs-handling` — `flow-cytometry/fcs-handling`
- `bio-flow-cytometry-gating-analysis` — `flow-cytometry/gating-analysis`

### gene-regulatory-networks

- `bio-gene-regulatory-networks-coexpression-networks` — `gene-regulatory-networks/coexpression-networks`
- `bio-gene-regulatory-networks-differential-networks` — `gene-regulatory-networks/differential-networks`
- `bio-gene-regulatory-networks-grn-inference` — `gene-regulatory-networks/grn-inference`
- `bio-gene-regulatory-networks-multiomics-grn` — `gene-regulatory-networks/multiomics-grn`
- `bio-gene-regulatory-networks-perturbation-simulation` — `gene-regulatory-networks/perturbation-simulation`
- `bio-gene-regulatory-networks-scenic-regulons` — `gene-regulatory-networks/scenic-regulons`

### genome-annotation

- `bio-genome-annotation-annotation-qc` — `genome-annotation/annotation-qc`
- `bio-genome-annotation-annotation-transfer` — `genome-annotation/annotation-transfer`
- `bio-genome-annotation-eukaryotic-gene-prediction` — `genome-annotation/eukaryotic-gene-prediction`
- `bio-genome-annotation-functional-annotation` — `genome-annotation/functional-annotation`
- `bio-genome-annotation-ncrna-annotation` — `genome-annotation/ncrna-annotation`
- `bio-genome-annotation-prokaryotic-annotation` — `genome-annotation/prokaryotic-annotation`
- `bio-genome-annotation-repeat-annotation` — `genome-annotation/repeat-annotation`

### genome-assembly

- `bio-genome-assembly-assembly-polishing` — `genome-assembly/assembly-polishing`
- `bio-genome-assembly-assembly-qc` — `genome-assembly/assembly-qc`
- `bio-genome-assembly-contamination-detection` — `genome-assembly/contamination-detection`
- `bio-genome-assembly-genome-profiling` — `genome-assembly/genome-profiling`
- `bio-genome-assembly-hifi-assembly` — `genome-assembly/hifi-assembly`
- `bio-genome-assembly-long-read-assembly` — `genome-assembly/long-read-assembly`
- `bio-genome-assembly-metagenome-assembly` — `genome-assembly/metagenome-assembly`
- `bio-genome-assembly-scaffolding` — `genome-assembly/scaffolding`
- `bio-genome-assembly-short-read-assembly` — `genome-assembly/short-read-assembly`

### genome-engineering

- `bio-genome-engineering-base-editing-design` — `genome-engineering/base-editing-design`
- `bio-genome-engineering-grna-design` — `genome-engineering/grna-design`
- `bio-genome-engineering-hdr-template-design` — `genome-engineering/hdr-template-design`
- `bio-genome-engineering-off-target-prediction` — `genome-engineering/off-target-prediction`
- `bio-genome-engineering-prime-editing-design` — `genome-engineering/prime-editing-design`

### genome-intervals

- `bio-genome-intervals-bed-file-basics` — `genome-intervals/bed-file-basics`
- `bio-genome-intervals-bedgraph-handling` — `genome-intervals/bedgraph-handling`
- `bio-genome-intervals-bigwig-tracks` — `genome-intervals/bigwig-tracks`
- `bio-genome-intervals-coverage-analysis` — `genome-intervals/coverage-analysis`
- `bio-genome-intervals-gtf-gff-handling` — `genome-intervals/gtf-gff-handling`
- `bio-genome-intervals-interval-arithmetic` — `genome-intervals/interval-arithmetic`
- `bio-genome-intervals-overlap-significance` — `genome-intervals/overlap-significance`
- `bio-genome-intervals-proximity-operations` — `genome-intervals/proximity-operations`

### hi-c-analysis

- `bio-hi-c-analysis-compartment-analysis` — `hi-c-analysis/compartment-analysis`
- `bio-hi-c-analysis-contact-pairs` — `hi-c-analysis/contact-pairs`
- `bio-hi-c-analysis-hic-data-io` — `hi-c-analysis/hic-data-io`
- `bio-hi-c-analysis-hic-differential` — `hi-c-analysis/hic-differential`
- `bio-hi-c-analysis-hic-visualization` — `hi-c-analysis/hic-visualization`
- `bio-hi-c-analysis-hichip-plac-loops` — `hi-c-analysis/hichip-plac-loops`
- `bio-hi-c-analysis-loop-calling` — `hi-c-analysis/loop-calling`
- `bio-hi-c-analysis-matrix-operations` — `hi-c-analysis/matrix-operations`
- `bio-hi-c-analysis-tad-detection` — `hi-c-analysis/tad-detection`

### imaging-mass-cytometry

- `bio-imaging-mass-cytometry-cell-segmentation` — `imaging-mass-cytometry/cell-segmentation`
- `bio-imaging-mass-cytometry-data-preprocessing` — `imaging-mass-cytometry/data-preprocessing`
- `bio-imaging-mass-cytometry-differential-analysis` — `imaging-mass-cytometry/differential-analysis`
- `bio-imaging-mass-cytometry-interactive-annotation` — `imaging-mass-cytometry/interactive-annotation`
- `bio-imaging-mass-cytometry-phenotyping` — `imaging-mass-cytometry/phenotyping`
- `bio-imaging-mass-cytometry-quality-metrics` — `imaging-mass-cytometry/quality-metrics`
- `bio-imaging-mass-cytometry-spatial-analysis` — `imaging-mass-cytometry/spatial-analysis`

### immunoinformatics

- `bio-immunoinformatics-epitope-prediction` — `immunoinformatics/epitope-prediction`
- `bio-immunoinformatics-immunogenicity-scoring` — `immunoinformatics/immunogenicity-scoring`
- `bio-immunoinformatics-mhc-binding-prediction` — `immunoinformatics/mhc-binding-prediction`
- `bio-immunoinformatics-mhc-class-ii-prediction` — `immunoinformatics/mhc-class-ii-prediction`
- `bio-immunoinformatics-neoantigen-prediction` — `immunoinformatics/neoantigen-prediction`
- `bio-immunoinformatics-tcr-epitope-binding` — `immunoinformatics/tcr-epitope-binding`

### liquid-biopsy

- `bio-analytical-validation` — `liquid-biopsy/analytical-validation`
- `bio-cfdna-preprocessing` — `liquid-biopsy/cfdna-preprocessing`
- `bio-ctdna-mutation-detection` — `liquid-biopsy/ctdna-mutation-detection`
- `bio-fragment-analysis` — `liquid-biopsy/fragment-analysis`
- `bio-longitudinal-monitoring` — `liquid-biopsy/longitudinal-monitoring`
- `bio-methylation-based-detection` — `liquid-biopsy/methylation-based-detection`
- `bio-tumor-fraction-estimation` — `liquid-biopsy/tumor-fraction-estimation`

### long-read-sequencing

- `bio-long-read-sequencing-basecalling` — `long-read-sequencing/basecalling`
- `bio-long-read-sequencing-clair3-variants` — `long-read-sequencing/clair3-variants`
- `bio-long-read-sequencing-haplotype-phasing` — `long-read-sequencing/haplotype-phasing`
- `bio-long-read-sequencing-isoseq-analysis` — `long-read-sequencing/isoseq-analysis`
- `bio-long-read-sequencing-long-read-alignment` — `long-read-sequencing/long-read-alignment`
- `bio-long-read-sequencing-long-read-qc` — `long-read-sequencing/long-read-qc`
- `bio-long-read-sequencing-medaka-polishing` — `long-read-sequencing/medaka-polishing`
- `bio-long-read-sequencing-nanopore-methylation` — `long-read-sequencing/nanopore-methylation`
- `bio-long-read-sequencing-structural-variants` — `long-read-sequencing/structural-variants`

### machine-learning

- `bio-machine-learning-atlas-mapping` — `machine-learning/atlas-mapping`
- `bio-machine-learning-biomarker-discovery` — `machine-learning/biomarker-discovery`
- `bio-machine-learning-omics-classifiers` — `machine-learning/omics-classifiers`
- `bio-machine-learning-survival-analysis` — `machine-learning/survival-analysis`

### metagenomics

- `bio-metagenomics-abundance` — `metagenomics/abundance-estimation`
- `bio-metagenomics-amr-detection` — `metagenomics/amr-detection`
- `bio-metagenomics-contamination-controls` — `metagenomics/contamination-controls`
- `bio-metagenomics-functional-profiling` — `metagenomics/functional-profiling`
- `bio-metagenomics-kraken` — `metagenomics/kraken-classification`
- `bio-metagenomics-metaphlan` — `metagenomics/metaphlan-profiling`
- `bio-metagenomics-strain-tracking` — `metagenomics/strain-tracking`
- `bio-metagenomics-visualization` — `metagenomics/metagenome-visualization`

### methylation-analysis

- `bio-methylation-array-preprocessing` — `methylation-analysis/array-preprocessing`
- `bio-methylation-array-qc-filtering` — `methylation-analysis/array-qc-filtering`
- `bio-methylation-bismark-alignment` — `methylation-analysis/bismark-alignment`
- `bio-methylation-calling` — `methylation-analysis/methylation-calling`
- `bio-methylation-cell-type-deconvolution` — `methylation-analysis/cell-type-deconvolution`
- `bio-methylation-differential-cpg` — `methylation-analysis/differential-cpg-testing`
- `bio-methylation-dmr-detection` — `methylation-analysis/dmr-detection`
- `bio-methylation-epigenetic-clocks` — `methylation-analysis/epigenetic-clocks`
- `bio-methylation-ewas-design` — `methylation-analysis/ewas-design`
- `bio-methylation-methylkit` — `methylation-analysis/methylkit-analysis`

### microbiome

- `bio-microbiome-functional-prediction` — `microbiome/functional-prediction`
- `bio-microbiome-qiime2-workflow` — `microbiome/qiime2-workflow`
- `bio-microbiome-taxonomy-assignment` — `microbiome/taxonomy-assignment`

### multi-omics-integration

- `bio-multi-omics-data-harmonization` — `multi-omics-integration/data-harmonization`
- `bio-multi-omics-integration-design` — `multi-omics-integration/integration-design`
- `bio-multi-omics-mixomics-analysis` — `multi-omics-integration/mixomics-analysis`
- `bio-multi-omics-mofa-integration` — `multi-omics-integration/mofa-integration`
- `bio-multi-omics-similarity-network` — `multi-omics-integration/similarity-network`

### phasing-imputation

- `bio-phasing-imputation-genotype-imputation` — `phasing-imputation/genotype-imputation`
- `bio-phasing-imputation-haplotype-phasing` — `phasing-imputation/haplotype-phasing`
- `bio-phasing-imputation-imputation-qc` — `phasing-imputation/imputation-qc`
- `bio-phasing-imputation-reference-panels` — `phasing-imputation/reference-panels`

### population-genetics

- `bio-population-genetics-association-testing` — `population-genetics/association-testing`
- `bio-population-genetics-linkage-disequilibrium` — `population-genetics/linkage-disequilibrium`
- `bio-population-genetics-plink-basics` — `population-genetics/plink-basics`
- `bio-population-genetics-population-structure` — `population-genetics/population-structure`
- `bio-population-genetics-scikit-allel-analysis` — `population-genetics/scikit-allel-analysis`
- `bio-population-genetics-selection-statistics` — `population-genetics/selection-statistics`

### primer-design

- `bio-primer-design-primer-basics` — `primer-design/primer-basics`
- `bio-primer-design-primer-specificity` — `primer-design/primer-specificity`
- `bio-primer-design-primer-validation` — `primer-design/primer-validation`
- `bio-primer-design-qpcr-primers` — `primer-design/qpcr-primers`

### proteomics

- `bio-proteomics-spectral-libraries` — `proteomics/spectral-libraries`

### read-alignment

- `bio-read-alignment-bowtie2-alignment` — `read-alignment/bowtie2-alignment`
- `bio-read-alignment-bwa-alignment` — `read-alignment/bwa-alignment`
- `bio-read-alignment-hisat2-alignment` — `read-alignment/hisat2-alignment`
- `bio-read-alignment-star-alignment` — `read-alignment/star-alignment`

### read-qc

- `bio-read-qc-adapter-trimming` — `read-qc/adapter-trimming`
- `bio-read-qc-contamination-screening` — `read-qc/contamination-screening`
- `bio-read-qc-fastp-workflow` — `read-qc/fastp-workflow`
- `bio-read-qc-quality-filtering` — `read-qc/quality-filtering`
- `bio-read-qc-quality-reports` — `read-qc/quality-reports`
- `bio-read-qc-rnaseq-qc` — `read-qc/rnaseq-qc`
- `bio-read-qc-umi-processing` — `read-qc/umi-processing`

### reporting

- `bio-reporting-automated-qc-reports` — `reporting/automated-qc-reports`
- `bio-reporting-figure-export` — `reporting/figure-export`
- `bio-reporting-jupyter-reports` — `reporting/jupyter-reports`
- `bio-reporting-publication-tables` — `reporting/publication-tables`
- `bio-reporting-quarto-reports` — `reporting/quarto-reports`
- `bio-reporting-rmarkdown-reports` — `reporting/rmarkdown-reports`

### restriction-analysis

- `bio-restriction-enzyme-selection` — `restriction-analysis/enzyme-selection`
- `bio-restriction-fragment-analysis` — `restriction-analysis/fragment-analysis`
- `bio-restriction-golden-gate-assembly` — `restriction-analysis/golden-gate-assembly`
- `bio-restriction-mapping` — `restriction-analysis/restriction-mapping`
- `bio-restriction-sites` — `restriction-analysis/restriction-sites`

### ribo-seq

- `bio-ribo-seq-initiation-site-mapping` — `ribo-seq/initiation-site-mapping`
- `bio-ribo-seq-orf-detection` — `ribo-seq/orf-detection`
- `bio-ribo-seq-riboseq-preprocessing` — `ribo-seq/riboseq-preprocessing`
- `bio-ribo-seq-ribosome-periodicity` — `ribo-seq/ribosome-periodicity`
- `bio-ribo-seq-ribosome-stalling` — `ribo-seq/ribosome-stalling`
- `bio-ribo-seq-translation-efficiency` — `ribo-seq/translation-efficiency`

### rna-quantification

- `bio-rna-quantification-alignment-free-quant` — `rna-quantification/alignment-free-quant`
- `bio-rna-quantification-count-matrix-qc` — `rna-quantification/count-matrix-qc`
- `bio-rna-quantification-featurecounts-counting` — `rna-quantification/featurecounts-counting`
- `bio-rna-quantification-tximport-workflow` — `rna-quantification/tximport-workflow`

### rna-structure

- `bio-rna-structure-covariation-analysis` — `rna-structure/covariation-analysis`
- `bio-rna-structure-ncrna-search` — `rna-structure/ncrna-search`
- `bio-rna-structure-secondary-structure-prediction` — `rna-structure/secondary-structure-prediction`
- `bio-rna-structure-structure-probing` — `rna-structure/structure-probing`

### sequence-io

- `bio-batch-processing` — `sequence-io/batch-processing`
- `bio-compressed-files` — `sequence-io/compressed-files`
- `bio-fastq-quality` — `sequence-io/fastq-quality`
- `bio-filter-sequences` — `sequence-io/filter-sequences`
- `bio-format-conversion` — `sequence-io/format-conversion`
- `bio-paired-end-fastq` — `sequence-io/paired-end-fastq`
- `bio-read-sequences` — `sequence-io/read-sequences`
- `bio-sequence-statistics` — `sequence-io/sequence-statistics`
- `bio-write-sequences` — `sequence-io/write-sequences`

### sequence-manipulation

- `bio-codon-usage` — `sequence-manipulation/codon-usage`
- `bio-motif-search` — `sequence-manipulation/motif-search`
- `bio-reverse-complement` — `sequence-manipulation/reverse-complement`
- `bio-seq-objects` — `sequence-manipulation/seq-objects`
- `bio-sequence-properties` — `sequence-manipulation/sequence-properties`
- `bio-sequence-slicing` — `sequence-manipulation/sequence-slicing`
- `bio-transcription-translation` — `sequence-manipulation/transcription-translation`

### single-cell

- `bio-single-cell-cell-communication` — `single-cell/cell-communication`
- `bio-single-cell-cnv-inference` — `single-cell/cnv-inference`
- `bio-single-cell-data-io` — `single-cell/data-io`
- `bio-single-cell-hashing-demultiplexing` — `single-cell/hashing-demultiplexing`
- `bio-single-cell-lineage-tracing` — `single-cell/lineage-tracing`
- `bio-single-cell-metabolite-communication` — `single-cell/metabolite-communication`
- `bio-single-cell-multimodal-integration` — `single-cell/multimodal-integration`
- `bio-single-cell-perturb-seq` — `single-cell/perturb-seq`
- `bio-single-cell-scatac-analysis` — `single-cell/scatac-analysis`
- `bio-single-cell-trajectory-inference` — `single-cell/trajectory-inference`

### small-rna-seq

- `bio-small-rna-seq-differential-mirna` — `small-rna-seq/differential-mirna`
- `bio-small-rna-seq-mirdeep2-analysis` — `small-rna-seq/mirdeep2-analysis`
- `bio-small-rna-seq-mirge3-analysis` — `small-rna-seq/mirge3-analysis`
- `bio-small-rna-seq-smrna-preprocessing` — `small-rna-seq/smrna-preprocessing`
- `bio-small-rna-seq-target-prediction` — `small-rna-seq/target-prediction`
- `bio-small-rna-seq-trf-pirna-profiling` — `small-rna-seq/trf-pirna-profiling`

### spatial-transcriptomics

- `bio-spatial-transcriptomics-high-resolution-binning` — `spatial-transcriptomics/high-resolution-binning`
- `bio-spatial-transcriptomics-image-analysis` — `spatial-transcriptomics/image-analysis`
- `bio-spatial-transcriptomics-spatial-communication` — `spatial-transcriptomics/spatial-communication`
- `bio-spatial-transcriptomics-spatial-data-io` — `spatial-transcriptomics/spatial-data-io`
- `bio-spatial-transcriptomics-spatial-deconvolution` — `spatial-transcriptomics/spatial-deconvolution`
- `bio-spatial-transcriptomics-spatial-domains` — `spatial-transcriptomics/spatial-domains`
- `bio-spatial-transcriptomics-spatial-multiomics` — `spatial-transcriptomics/spatial-multiomics`
- `bio-spatial-transcriptomics-spatial-neighbors` — `spatial-transcriptomics/spatial-neighbors`
- `bio-spatial-transcriptomics-spatial-preprocessing` — `spatial-transcriptomics/spatial-preprocessing`
- `bio-spatial-transcriptomics-spatial-proteomics` — `spatial-transcriptomics/spatial-proteomics`
- `bio-spatial-transcriptomics-spatial-statistics` — `spatial-transcriptomics/spatial-statistics`
- `bio-spatial-transcriptomics-spatial-visualization` — `spatial-transcriptomics/spatial-visualization`

### structural-biology

- `bio-structural-biology-alphafold-predictions` — `structural-biology/alphafold-predictions`
- `bio-structural-biology-binding-site-detection` — `structural-biology/binding-site-detection`
- `bio-structural-biology-geometric-analysis` — `structural-biology/geometric-analysis`
- `bio-structural-biology-interface-analysis` — `structural-biology/interface-analysis`
- `bio-structural-biology-modern-structure-prediction` — `structural-biology/modern-structure-prediction`
- `bio-structural-biology-structure-io` — `structural-biology/structure-io`
- `bio-structural-biology-structure-modification` — `structural-biology/structure-modification`
- `bio-structural-biology-structure-navigation` — `structural-biology/structure-navigation`
- `bio-structural-biology-structure-preparation` — `structural-biology/structure-preparation`
- `bio-structural-biology-structure-validation` — `structural-biology/structure-validation`

### systems-biology

- `bio-systems-biology-community-metabolic-modeling` — `systems-biology/community-metabolic-modeling`
- `bio-systems-biology-context-specific-models` — `systems-biology/context-specific-models`
- `bio-systems-biology-flux-balance-analysis` — `systems-biology/flux-balance-analysis`
- `bio-systems-biology-gene-essentiality` — `systems-biology/gene-essentiality`
- `bio-systems-biology-metabolic-reconstruction` — `systems-biology/metabolic-reconstruction`
- `bio-systems-biology-model-curation` — `systems-biology/model-curation`
- `bio-systems-biology-strain-design` — `systems-biology/strain-design`

### tcr-bcr-analysis

- `bio-tcr-bcr-analysis-immcantation-analysis` — `tcr-bcr-analysis/immcantation-analysis`
- `bio-tcr-bcr-analysis-mixcr-analysis` — `tcr-bcr-analysis/mixcr-analysis`
- `bio-tcr-bcr-analysis-repertoire-visualization` — `tcr-bcr-analysis/repertoire-visualization`
- `bio-tcr-bcr-analysis-scirpy-analysis` — `tcr-bcr-analysis/scirpy-analysis`
- `bio-tcr-bcr-analysis-specificity-annotation` — `tcr-bcr-analysis/specificity-annotation`
- `bio-tcr-bcr-analysis-vdjtools-analysis` — `tcr-bcr-analysis/vdjtools-analysis`

### temporal-genomics

- `bio-temporal-genomics-circadian-rhythms` — `temporal-genomics/circadian-rhythms`
- `bio-temporal-genomics-differential-rhythmicity` — `temporal-genomics/differential-rhythmicity`
- `bio-temporal-genomics-periodicity-detection` — `temporal-genomics/periodicity-detection`
- `bio-temporal-genomics-temporal-clustering` — `temporal-genomics/temporal-clustering`
- `bio-temporal-genomics-temporal-grn` — `temporal-genomics/temporal-grn`
- `bio-temporal-genomics-trajectory-modeling` — `temporal-genomics/trajectory-modeling`

### variant-calling

- `bio-consensus-sequences` — `variant-calling/consensus-sequences`
- `bio-gatk-variant-calling` — `variant-calling/gatk-variant-calling`
- `bio-variant-calling` — `variant-calling/variant-calling`
- `bio-variant-calling-clinical-interpretation` — `variant-calling/clinical-interpretation`
- `bio-variant-calling-deepvariant` — `variant-calling/deepvariant`
- `bio-variant-calling-joint-calling` — `variant-calling/joint-calling`
- `bio-variant-calling-structural-variant-calling` — `variant-calling/structural-variant-calling`

### workflow-management

- `bio-workflow-management-cwl-workflows` — `workflow-management/cwl-workflows`
- `bio-workflow-management-nextflow-pipelines` — `workflow-management/nextflow-pipelines`
- `bio-workflow-management-nf-core-pipelines` — `workflow-management/nf-core-pipelines`
- `bio-workflow-management-snakemake-workflows` — `workflow-management/snakemake-workflows`
- `bio-workflow-management-wdl-workflows` — `workflow-management/wdl-workflows`

### workflows

- `bio-workflows-atacseq-pipeline` — `workflows/atacseq-pipeline`
- `bio-workflows-biomarker-pipeline` — `workflows/biomarker-pipeline`
- `bio-workflows-causal-genomics-pipeline` — `workflows/causal-genomics-pipeline`
- `bio-workflows-chipseq-pipeline` — `workflows/chipseq-pipeline`
- `bio-workflows-clinical-trial-pipeline` — `workflows/clinical-trial-pipeline`
- `bio-workflows-clip-pipeline` — `workflows/clip-pipeline`
- `bio-workflows-cnv-pipeline` — `workflows/cnv-pipeline`
- `bio-workflows-crispr-editing-pipeline` — `workflows/crispr-editing-pipeline`
- `bio-workflows-cytometry-pipeline` — `workflows/cytometry-pipeline`
- `bio-workflows-edna-pipeline` — `workflows/edna-pipeline`
- `bio-workflows-expression-to-pathways` — `workflows/expression-to-pathways`
- `bio-workflows-fastq-to-variants` — `workflows/fastq-to-variants`
- `bio-workflows-genome-annotation-pipeline` — `workflows/genome-annotation-pipeline`
- `bio-workflows-genome-assembly-pipeline` — `workflows/genome-assembly-pipeline`
- `bio-workflows-grn-pipeline` — `workflows/grn-pipeline`
- `bio-workflows-gwas-pipeline` — `workflows/gwas-pipeline`
- `bio-workflows-hic-pipeline` — `workflows/hic-pipeline`
- `bio-workflows-imc-pipeline` — `workflows/imc-pipeline`
- `bio-workflows-liquid-biopsy-pipeline` — `workflows/liquid-biopsy-pipeline`
- `bio-workflows-longread-sv-pipeline` — `workflows/longread-sv-pipeline`
- `bio-workflows-merip-pipeline` — `workflows/merip-pipeline`
- `bio-workflows-metabolic-modeling-pipeline` — `workflows/metabolic-modeling-pipeline`
- `bio-workflows-metagenomics-pipeline` — `workflows/metagenomics-pipeline`
- `bio-workflows-methylation-pipeline` — `workflows/methylation-pipeline`
- `bio-workflows-microbiome-pipeline` — `workflows/microbiome-pipeline`
- `bio-workflows-multi-omics-pipeline` — `workflows/multi-omics-pipeline`
- `bio-workflows-multiome-pipeline` — `workflows/multiome-pipeline`
- `bio-workflows-neoantigen-pipeline` — `workflows/neoantigen-pipeline`
- `bio-workflows-outbreak-pipeline` — `workflows/outbreak-pipeline`
- `bio-workflows-riboseq-pipeline` — `workflows/riboseq-pipeline`
- `bio-workflows-rnaseq-to-de` — `workflows/rnaseq-to-de`
- `bio-workflows-smrna-pipeline` — `workflows/smrna-pipeline`
- `bio-workflows-somatic-variant-pipeline` — `workflows/somatic-variant-pipeline`
- `bio-workflows-spatial-pipeline` — `workflows/spatial-pipeline`
- `bio-workflows-splicing-pipeline` — `workflows/splicing-pipeline`
- `bio-workflows-tcr-pipeline` — `workflows/tcr-pipeline`
- `bio-workflows-timecourse-pipeline` — `workflows/timecourse-pipeline`
