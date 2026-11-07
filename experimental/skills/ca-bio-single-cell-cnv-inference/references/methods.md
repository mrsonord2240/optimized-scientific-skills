# Source methods and limitations

Source methods and package APIs are unverified. Do not infer that SCEVAN classification succeeded merely because an error was caught; inspect a successfully returned result. An arbitrary median split is not a validated malignant classifier.

## Governing principle

Averaged expression over a genomic window is a PROXY for DNA copy number, not a measurement of it. A chromosome-arm gain raises the average expression of the many genes sitting on that arm, and a loss lowers it, so smoothing expression across long runs of contiguous genes reconstructs a coarse copy-number landscape. Three consequences follow and drive every decision. First, the signal is noisy and indirect: it must be smoothed over windows of dozens to hundreds of neighboring genes, and its resolution is chromosome-arm or large-segment (around 5 Mb), never focal genes or exons - this is the hard line separating it from DNA-based copy-number (copy-number/cnvkit-analysis, copy-number/gatk-cnv), which measures DNA read depth and resolves focal amplifications and deletions. Second, the value is RELATIVE: copy-number is only defined against a copy-neutral baseline, so a NORMAL reference is required, and the choice of that reference is the single most consequential decision - a wrong or mismatched reference fabricates CNVs out of ordinary cell-type expression differences. Third, calling a cell malignant is a CLUSTERING decision on this noisy proxy, so it is a hypothesis that needs orthogonal support (lineage markers, mutations, allele evidence), and subclone calls are even softer hypotheses.

Tumor CNV profiles are patient-PRIVATE, so CNV inference runs PER SAMPLE / per patient. Integrating or batch-correcting cells across patients before inference mixes distinct private karyotypes and erases the per-patient signal the analysis depends on (single-cell/batch-integration notes this). Integrate across patients only for shared transcriptional-state analysis, never as input to CNV calling.

CNV-quiet tumors are invisible to this approach. Many hematologic and low-grade tumors carry little large-scale CNV, so expression-based inference returns a near-flat profile - absence of inferred CNV is NOT evidence that cells are normal, and a confident malignant call then requires allele evidence or orthogonal markers. Balanced whole-genome doubling is invisible for the same reason and to all expression methods, not just copyKAT: per-cell library-size normalization removes a uniform ploidy multiple, so a 4N WGD reads copy-neutral against a 2N reference for inferCNV and Numbat's expression channel alike, and only allele or ploidy evidence reveals it.

## Choosing a CNV-inference method

| Method | Model / assumption | Reference | Allele-aware | Use when | Fails when |
|--------|--------------------|-----------|--------------|----------|------------|
| inferCNV | Smooth expression along chromosomes vs a normal reference; optional HMM for discrete states | Reference-based (needs normal cells) | No | A clean in-sample normal reference exists; want interpretable arm-level heatmap + HMM states | No trustworthy reference; CNV-quiet tumor; needs focal resolution |
| copyKAT | Bayesian segmentation + hierarchical clustering + GMM; estimates diploid baseline internally | Reference-free (optional known normals) | No | No reference cells annotated; want a quick aneuploid-vs-diploid call at ~5 Mb | Mostly-aneuploid sample with no diploid baseline; whole-genome doubling confuses the root |
| SCEVAN | Variational multichannel segmentation sharing breakpoints across a clone | Reference-free | No | Want automatic malignant/non-malignant + subclones in one call | Same baseline ambiguity as copyKAT; very sparse data |
| Numbat | Joint expression + phased B-allele frequency + population haplotypes; iterative phylogeny | Reference (expression) + population phasing | Yes | Best subclone resolution needed; copy-neutral LOH matters; allele counts obtainable | No BAM/phasing available; very low SNP coverage (shallow or snRNA) |

Reference-based (inferCNV, Numbat expression side) is the most reliable when a clean normal reference is available; reference-free (copyKAT, SCEVAN) trades that for not needing one but is vulnerable to baseline ambiguity. Expression-only methods (inferCNV, copyKAT, SCEVAN) are simpler and need only counts; the allele-aware method (Numbat) is the most powerful for subclones and copy-neutral events but requires per-cell allele counts and phasing. Run an expression method to get the malignant/normal split, then Numbat when subclone structure is the question, and reconcile. When methods compete, verify current best practice against the installed package documentation before committing to one.

## Choosing the normal reference

**Goal:** Pick a copy-neutral baseline that defines what "no CNV" looks like, since the reference choice determines whether the inferred CNVs are real or artifacts.

**Approach:** Prefer non-malignant cells from the SAME sample (T cells, B cells, myeloid, endothelial, fibroblasts identified by lineage markers) because they share the patient, protocol, and ambient-RNA background; fall back to an external normal only when no in-sample normals exist, and expect batch artifacts.

The reference cells must be confidently non-malignant and abundant enough for a stable baseline mean: aim for tens to hundreds of reference cells, not a handful, because a few cells give a noisy baseline that fabricates CNVs in every observation cell. A reference that is itself a malignant or stressed population, or a single mismatched cell type, will make every other cell look aneuploid relative to it. The reference must also MATCH the tumor's sex: because inferCNV works on expression, X-inactivation dosage-compensates most chrX expression so there is no uniform two-fold chrX shift, but an external or cross-individual normal (including shipped or pooled-donor references) of the opposite sex still fabricates a convincing uniform sex-chromosome CNV through chrY genes present in XY versus near-absent in XX, XIST high in XX versus silent in XY, and the attenuated X-inactivation escape genes - match sexes or drop chrX and chrY before inference. With no annotated normals, use a reference-free method (copyKAT, SCEVAN) rather than guessing a reference; do not pass tumor-contaminated cells as the reference.

## Numbat - haplotype-aware allele + expression

**Goal:** Resolve subclones and copy-neutral LOH by combining smoothed expression with phased B-allele frequencies.

**Approach:** Generate per-cell allele counts and population phasing with the `pileup_and_phase.R` preprocessing script, build an expression reference from matched normals (or the shipped `ref_hca`), then run the joint model.

```bash
Rscript pileup_and_phase.R --label tumor1 --samples tumor1 \
    --bams tumor1.bam --barcodes barcodes.tsv \
    --gmap genetic_map_hg38_withX.txt.gz --snpvcf genome1K.phase3.SNP.vcf \
    --paneldir 1000G_panel/ --outdir numbat_out/ --ncores 8
```

```r
library(numbat)

ref <- aggregate_counts(count_mat_normal, cell_annot)
out <- run_numbat(
    count_mat,
    lambdas_ref = ref,
    df_allele = df_allele,
    genome = 'hg38',
    t = 1e-5,
    ncores = 4,
    plot = TRUE,
    out_dir = 'numbat_out')

nb <- Numbat$new(out_dir = 'numbat_out')
```

`df_allele` is the allele dataframe written by `pileup_and_phase.R`, with one row per SNP per cell. `run_numbat()` hard-requires ten input columns: `cell`, `snp_id`, `CHROM`, `POS`, `cM`, `REF`, `ALT`, `AD`, `DP`, `GT` - a `df_allele` missing `cM`, `REF`, or `ALT` is rejected before any computation runs. `check_allele_df()` itself lists an eleventh, `gene`, but `run_numbat()` calls `annotate_genes(df_allele, gtf)` first, which adds `gene` from the SNP position and the genome's gene table (and overwrites any `gene` you pass), so do not build it by hand; only a direct call to `numbat:::check_allele_df()` on an unannotated frame fails on `gene` (checked on numbat 1.5.2: a ten-column frame passes `check_allele_df(annotate_genes(df, gtf_hg38))`, and fails on `gene` without the annotate step). `cM` (genetic-map position) comes from the genetic-map file passed to `pileup_and_phase.R` via `--gmap`; `REF`/`ALT` are the reference/alternate alleles at each SNP. `lambdas_ref` is a gene-by-cell-type expression reference from `aggregate_counts(count_mat, cell_annot)` where `cell_annot` has `cell` and `group` columns, or the package-shipped `ref_hca`. `t` is the HMM transition probability. The loaded `Numbat` object exposes `clone_post` (clone assignments) and per-cell copy-number posteriors. Numbat needs no paired-normal DNA but does need a BAM and phasing reference.

## Threshold and parameter reference

| Parameter | Tool | Default / typical | Rationale |
|-----------|------|-------------------|-----------|
| cutoff | inferCNV | 0.1 (droplet), 1 (Smart-seq) | Drops genes below mean expression; droplet data are sparse so the threshold is lower |
| HMM_type | inferCNV | i6 | Six copy states (0 to 3+) give finer state calls; i3 (del/neutral/amp) is simpler and more robust |
| ref_group_names | inferCNV | named normals | The copy-neutral baseline; NULL only when no reference exists, at a reliability cost |
| ngene.chr | copyKAT | 5 | Minimum genes per chromosome per cell so a cell has signal on each chromosome |
| win.size | copyKAT | 25 | Genes per segment for smoothing; larger windows denoise but blur small events |
| KS.cut | copyKAT | 0.1 | Segmentation stringency; raise for fewer, larger segments in noisy data |
| genome | copyKAT | hg20 | Must match the assembly the gene coordinates came from (hg20 or mm10) |
| min.gene.per.cell | copyKAT | 200 | Minimum genes detected per cell to keep it; lower for small/custom gene panels (default drops every cell on a panel under ~200 genes) |
| t | Numbat | 1e-5 | HMM transition probability; lower favors longer segments |
| df_allele columns | Numbat | required | `run_numbat()` requires `cell, snp_id, CHROM, POS, cM, REF, ALT, AD, DP, GT`; `gene` is added internally by `annotate_genes()` (see the Numbat section); a hand-built `df_allele` missing `cM`/`REF`/`ALT` is rejected before any computation runs |
| ngenes_chr | SCEVAN | 5 | Minimum genes per chromosome per cell, same role as copyKAT's `ngene.chr` |
| resolution (~5 Mb) | all expression methods | ~5 Mb | The floor of expression-based CNV; focal events below this are invisible |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Every cell looks aneuploid, including the immune cells | Reference cells are tumor-contaminated or a mismatched cell type | Re-pick confident non-malignant in-sample normals by lineage markers; do not use a stressed/malignant population as reference |
| Flat profile, no CNVs called, but cells are clearly tumor | CNV-quiet / low-grade tumor, or snRNA sparsity | Absence of CNV is not normality; switch to allele-aware Numbat or confirm malignancy with markers/mutations |
| copyKAT cannot find a diploid baseline / labels almost all aneuploid | Mostly-aneuploid sample with no internal diploid cells, or balanced whole-genome doubling that library-size normalization hides | Supply known-normal barcodes via `norm.cell.names`, or use inferCNV with an external reference; confirm ploidy with allele/DNA evidence since balanced WGD looks copy-neutral to all expression methods |
| Uniform chrX gain/loss and a chrY call across all tumor cells | Reference and tumor are opposite sex (external/shipped/pooled-donor reference) | Match reference and tumor sex, or drop chrX and chrY before inference |
| Recurrent localized segment over 6p, 14q, 22q, or 2p that tracks lineage | High, variable HLA (MHC 6p) and immunoglobulin/TCR expression clusters genomically and mimics a segment, especially with immune references | Mask or distrust segments over HLA and Ig/TCR loci; do not call them CNV |
| Stripe artifacts that track proliferating cells | Cell-cycle and high-expression gene programs mimic CNV | Account for / regress cell cycle, enable denoising, and treat cycle-correlated bands skeptically |
| Subclones from inferCNV/copyKAT do not replicate | Expression-only subclone calls are weak hypotheses | Confirm with Numbat allele evidence or DNA; report subclones as hypotheses |
| CNV signal vanishes after integrating samples | Cross-patient integration erased patient-private karyotypes | Run CNV inference per sample BEFORE any cross-patient integration |
| Genes silently dropped / wrong chromosome bands | Gene-order file or genome build mismatched to the counts | Match `gene_order_file` and `genome` to the same assembly and gene IDs as the matrix |
| SCEVAN `pipelineCNA()` errors after printing "found N tumor cells" | Its own internal plotting step (no `plot=FALSE` switch) fails on a sparse or small cohort, after classification already completed | Surface the error and inspect saved artifacts; a console count is not a recovered per-cell classification |
| SCEVAN "all cells are filtered" at step 5 | Gene panel doesn't span enough of the 22 autosomes SCEVAN checks per cell (`ngenes_chr` genes on each) | Use a genome-wide (or broadly multi-chromosome) panel, not a handful of loci; this is an algorithmic requirement, not tunable away with `ngenes_chr` alone |

## References

- Patel AP et al. 2014, Science 344:1396-1401 - single-cell RNA-seq of glioblastoma; first per-cell CNV estimation from averaged expression.
- Tirosh I et al. 2016, Science 352:189-196 - melanoma scRNA-seq; inferring large-scale CNV from smoothed expression to separate malignant cells (the inferCNV approach).
- Gao R et al. 2021, Nat Biotechnol 39:599-608 - copyKAT; reference-free Bayesian segmentation calling aneuploid vs diploid and subclones at ~5 Mb.
- Gao T et al. 2023, Nat Biotechnol 41:417-426 - Numbat; haplotype-aware joint allele and expression somatic CNV inference from scRNA-seq.
- Muller S et al. 2018, Bioinformatics 34:3217-3219 - CONICS/CONICSmat; arm-level CNV from expression mapped to tumor sub-clones.
- De Falco A et al. 2023, Nat Commun 14:1074 - SCEVAN; variational multichannel segmentation auto-classifying malignant cells and clonal substructure.
