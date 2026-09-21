#!/usr/bin/env Rscript
# Reference: FRASER 2.2.0 (Bioc 3.20, R 4.4) and 2.6.1 (Bioc 3.22, R 4.5) | Verify API if version differs
# FRASER 2 outlier splicing detection for rare-disease RESEARCH (not a clinical report).
# Usage: Rscript fraser2_rare_disease.R [bam_dir=bams] [patient_id=PATIENT_001] [working_dir=fraser_workdir]
#   bam_dir holds one coordinate-sorted, indexed <sampleID>.bam per sample (patient + controls, same tissue)
#
# Workflow: count split reads -> PSI values -> filter -> pick q -> fit (PCA) -> results for the patient

suppressPackageStartupMessages({
    library(FRASER)
    library(BiocParallel)
})

args <- commandArgs(TRUE)
bam_dir <- if (length(args) >= 1) args[1] else 'bams'
patient_id <- if (length(args) >= 2) args[2] else 'PATIENT_001'
working_dir <- if (length(args) >= 3) args[3] else 'fraser_workdir'
bp <- if (.Platform$OS.type == 'windows') SerialParam() else MulticoreParam(8)

bam_files <- list.files(bam_dir, pattern = '\\.bam$', full.names = TRUE)
sample_table <- data.frame(
    sampleID = sub('\\.bam$', '', basename(bam_files)),
    bamFile = bam_files,
    pairedEnd = TRUE
)
stopifnot(patient_id %in% sample_table$sampleID)
if (nrow(sample_table) < 20) {
    warning('Only ', nrow(sample_table), ' samples: FRASER missed a planted skipping event at n=12 and n=8 (see SKILL.md, Cohort Size)')
}

# colData must be an S4Vectors DataFrame; a plain data.frame errors on FRASER 2.2.0 and 2.6.1
fds <- FraserDataSet(
    colData = S4Vectors::DataFrame(sample_table),
    workingDir = working_dir,
    name = 'rare_disease_cohort'
)

fds <- countRNAData(fds, BPPARAM = bp)
fds <- calculatePSIValues(fds)

# FRASER defaults for the variability filter (quantile = 0.75, quantileMinExpression = 10)
fds <- filterExpressionAndVariability(fds, minExpressionInOneSample = 20, minDeltaPsi = 0.0)
cat(sprintf('%d junctions kept after filtering, %d samples\n', nrow(fds), ncol(fds)))

# FRASER 2 fits the Intron Jaccard Index only (the default metric)
fitMetrics(fds) <- 'jaccard'
currentType(fds) <- 'jaccard'

# q = size of the PCA latent space. Estimate it per cohort: estimateBestQ() is in FRASER 2.6.1 (optimal
# hard threshold, no injection); FRASER 2.2.0 only has optimHyperParams() (injected-outlier grid search).
est_q <- if (exists('estimateBestQ', envir = asNamespace('FRASER'), inherits = FALSE)) {
    FRASER::estimateBestQ
} else {
    FRASER::optimHyperParams
}
fds <- est_q(fds, type = 'jaccard', plot = FALSE)
q_jaccard <- bestQ(fds, type = 'jaccard')
cat(sprintf('q = %d for %d samples\n', q_jaccard, ncol(fds)))

# implementation = 'PCA' is the default and is deterministic; 'AE' is reproducible only with set.seed(1) AND BPPARAM = SerialParam(RNGseed = 1) together
fds <- FRASER(fds, q = c(jaccard = q_jaccard), implementation = 'PCA', BPPARAM = bp)
saveFraserDataSet(fds, dir = working_dir)

all_results <- results(
    fds,
    psiType = 'jaccard',
    padjCutoff = 0.05,
    deltaPsiCutoff = 0.1
)

patient_results <- as.data.frame(all_results)
if (nrow(patient_results) > 0) {    # an empty result has no columns
    patient_results <- patient_results[patient_results$sampleID == patient_id, ]
    patient_results <- patient_results[order(patient_results$padjust), ]
}

write.table(
    patient_results,
    file = file.path(working_dir, paste0(patient_id, '_outliers.tsv')),
    sep = '\t', quote = FALSE, row.names = FALSE
)

# Volcano plot for the patient (Rscript would otherwise write Rplots.pdf)
pdf(file.path(working_dir, paste0(patient_id, '_volcano.pdf')))
plotVolcano(fds, sampleID = patient_id, type = 'jaccard')
invisible(dev.off())

cat(sprintf('%d aberrant junctions in %s (padj<0.05, |delta|>=0.1); %d in all samples\n',
            nrow(patient_results), patient_id, length(all_results)))
