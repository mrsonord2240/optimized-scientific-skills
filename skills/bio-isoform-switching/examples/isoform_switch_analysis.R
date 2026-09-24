#!/usr/bin/env Rscript
# Isoform switch analysis with IsoformSwitchAnalyzeR: Salmon quant -> switches -> consequences -> plot.
# Checked on IsoformSwitchAnalyzeR 2.6.0 / Bioconductor 3.20, Salmon 2.7.0 output. Verify API if versions differ.
#
#   Rscript isoform_switch_analysis.R
#       Demo: builds a small synthetic Salmon dataset (planted switches) in a temp dir, runs everything,
#       and asserts that the planted switches are recovered.
#   Rscript isoform_switch_analysis.R <salmon_dir> <annotation.gtf> <transcripts.fa> <sample_metadata.tsv> [out_dir] [annotator_dir]
#       Real data. salmon_dir holds one sub-directory per sample, each with quant.sf.
#       sample_metadata.tsv: tab-separated, columns sample_id (= sub-directory name), condition, optional
#       covariate columns (batch, ...). Rows are joined to the quantification BY NAME, never by position.
#       annotator_dir (optional): cpc2_result.txt, pfam_scanfmt.txt (see hmmscan_to_pfamscan.py),
#       signalp_results.txt, iupred2a_result.txt; each present file adds its consequence type.

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))

# ---- toy data (demo mode only) ----------------------------------------------------------------
# 60 genes x 3 isoforms (A canonical; B = A + exon P; C = A minus exon E3), 3 v 3 samples.
# Genes 1-8: B rises 0.10 -> 0.55 and B carries a premature stop in P (NMD). Genes 9-14: C rises 0.20 -> 0.60.
# Sample IDs are deliberately not in condition order, so a positional design vector would mislabel them.
make_toy_data <- function(dir) {
    set.seed(11)
    nt <- c('A', 'C', 'G', 'T')
    sense <- setdiff(apply(expand.grid(nt, nt, nt), 1, paste, collapse = ''), c('TAA', 'TAG', 'TGA', 'ATG'))
    codons <- function(k) paste(sample(sense, k, replace = TRUE), collapse = '')
    utr <- function(n) repeat { s <- paste(sample(nt, n, replace = TRUE), collapse = ''); if (!grepl('ATG', s)) return(s) }
    G <- 60; gid <- sprintf('TOY%03d', 1:G)
    type <- c(rep('poison', 8), rep('skip', 6), rep('null', G - 14))
    strand <- sample(c('+', '-'), G, replace = TRUE)
    tx_seq <- character(); gtf <- character(); pos <- 1000L
    for (i in 1:G) {
        ex <- list(E1 = paste0(utr(30), 'ATG', codons(60)), E2 = codons(50),
                   P = if (type[i] == 'poison') paste0('TAA', codons(19)) else codons(20),
                   E3 = codons(40), E4 = paste0(codons(40), 'TAA', utr(150)))
        cur <- pos; coords <- list()
        for (e in names(ex)) { L <- nchar(ex[[e]]); coords[[e]] <- c(cur, cur + L - 1L); cur <- cur + L + sample(400:900, 1) }
        pos <- cur + 2000L
        if (strand[i] == '-') { tot <- max(unlist(coords)) + min(unlist(coords)); coords <- lapply(coords, function(cc) sort(tot - cc)) }
        defs <- list(A = c('E1', 'E2', 'E3', 'E4'), B = c('E1', 'E2', 'P', 'E3', 'E4'), C = c('E1', 'E2', 'E4'))
        for (t in names(defs)) {
            tid <- paste0(gid[i], '_', t); tx_seq[tid] <- paste(unlist(ex[defs[[t]]]), collapse = '')
            for (e in defs[[t]]) gtf <- c(gtf, sprintf('chr1\tsynth\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; gene_name "%s";',
                                                       coords[[e]][1], coords[[e]][2], strand[i], gid[i], tid, gid[i]))
        }
    }
    writeLines(gtf, file.path(dir, 'annotation.gtf'))
    writeLines(unlist(lapply(names(tx_seq), function(t) c(paste0('>', t), tx_seq[[t]]))), file.path(dir, 'transcripts.fa'))
    base <- matrix(c(0.60, 0.10, 0.30), G, 3, byrow = TRUE); alt <- base
    alt[type == 'poison', ] <- rep(c(0.25, 0.55, 0.20), each = sum(type == 'poison'))
    alt[type == 'skip', ]   <- rep(c(0.20, 0.20, 0.60), each = sum(type == 'skip'))
    meta <- data.frame(sample_id = c('SRR9001', 'SRR9002', 'SRR9003', 'SRR9004', 'SRR9005', 'SRR9006'),
                       condition = c('treatment', 'control', 'treatment', 'control', 'control', 'treatment'))
    write.table(meta, file.path(dir, 'sample_metadata.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    tlen <- nchar(tx_seq); efflen <- pmax(tlen - 74, 20)
    gene_mu <- exp(rnorm(G, 6.7, 0.4)); rdir <- function(a) { g <- rgamma(3, a); g / sum(g) }
    for (j in seq_len(nrow(meta))) {
        d <- file.path(dir, 'salmon', meta$sample_id[j]); dir.create(d, recursive = TRUE)
        n <- unlist(lapply(1:G, function(i) {
            p <- rdir(60 * (if (meta$condition[j] == 'treatment') alt[i, ] else base[i, ]))
            as.vector(rmultinom(1, rnbinom(1, mu = gene_mu[i], size = 30), p)) }))
        tpm <- (n / efflen) / sum(n / efflen) * 1e6
        write.table(data.frame(Name = names(tx_seq), Length = tlen, EffectiveLength = efflen, TPM = tpm, NumReads = n),
                    file.path(d, 'quant.sf'), sep = '\t', quote = FALSE, row.names = FALSE)
    }
    list(salmon_dir = file.path(dir, 'salmon'), gtf = file.path(dir, 'annotation.gtf'), fasta = file.path(dir, 'transcripts.fa'),
         metadata = file.path(dir, 'sample_metadata.tsv'), planted = gid[type != 'null'], poison = gid[type == 'poison'])
}

# ---- arguments ----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
demo <- length(args) == 0
if (demo) {
    work <- file.path(tempdir(), 'isoswitch_demo'); dir.create(work, recursive = TRUE)
    toy <- make_toy_data(work)
    salmon_dir <- toy$salmon_dir; gtf <- toy$gtf; fasta <- toy$fasta; meta_file <- toy$metadata
    out_dir <- file.path(work, 'out'); annot_dir <- NULL
} else {
    if (length(args) < 4) stop('usage: isoform_switch_analysis.R <salmon_dir> <annotation.gtf> <transcripts.fa> <sample_metadata.tsv> [out_dir] [annotator_dir]')
    salmon_dir <- args[1]; gtf <- args[2]; fasta <- args[3]; meta_file <- args[4]
    out_dir <- if (length(args) >= 5) args[5] else 'isoswitch_out'
    annot_dir <- if (length(args) >= 6) args[6] else NULL
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. import (raw NumReads as counts; see SKILL.md "Count route") ------------------------------
salmonQuant <- importIsoformExpression(parentDir = salmon_dir, calculateCountsFromAbundance = FALSE,
                                       addIsofomIdAsColumn = TRUE, showProgress = FALSE)

# Design joined to the quantification BY SAMPLE NAME. importIsoformExpression sorts samples alphabetically,
# so a hand-typed condition vector silently mislabels samples unless the IDs happen to sort into condition order.
meta <- read.delim(meta_file, stringsAsFactors = FALSE)
stopifnot('metadata needs columns sample_id and condition' = all(c('sample_id', 'condition') %in% colnames(meta)),
          'duplicated sample_id in metadata' = !anyDuplicated(meta$sample_id))
ids <- setdiff(colnames(salmonQuant$counts), 'isoform_id')
if (!setequal(ids, meta$sample_id))
    stop('sample IDs differ between quantification and metadata.\n  only in quant: ', paste(setdiff(ids, meta$sample_id), collapse = ', '),
         '\n  only in metadata: ', paste(setdiff(meta$sample_id, ids), collapse = ', '))
design <- meta[match(ids, meta$sample_id), , drop = FALSE]
colnames(design)[colnames(design) == 'sample_id'] <- 'sampleID'
design <- design[, c('sampleID', 'condition', setdiff(colnames(design), c('sampleID', 'condition'))), drop = FALSE]
rownames(design) <- NULL
stopifnot(identical(design$sampleID, ids))
nrep <- table(design$condition)
if (length(nrep) != 2 || min(nrep) < 2) stop('need exactly 2 conditions with >= 2 replicates each; got: ', paste(names(nrep), nrep, collapse = ', '))
if (min(nrep) < 3) warning('fewer than 3 replicates per condition: raw-count calls are not reliable (label-shuffled splits of a real 2 v 2 set were called almost as often as the true one); see SKILL.md "Null check"')
print(design)

# importRdata(detectUnwantedEffects = TRUE) may run stochastic sva. Keep the
# default enabled, but make a run reproducible and report this seed with results.
set.seed(1)
aSwitchList <- importRdata(isoformCountMatrix = salmonQuant$counts, isoformRepExpression = salmonQuant$abundance,
                           designMatrix = design, isoformExonAnnoation = gtf, isoformNtFasta = fasta,
                           addAnnotatedORFs = TRUE, showProgress = FALSE)

# ---- 2. filter and test (DEXSeq up to 5 replicates per condition, satuRn above) ------------------
aSwitchList <- preFilter(aSwitchList, geneExpressionCutoff = 1, isoformExpressionCutoff = 0, IFcutoff = 0.01,
                         removeSingleIsoformGenes = TRUE, keepIsoformInAllConditions = TRUE)
aSwitchList <- if (max(nrep) > 5) {
    isoformSwitchTestSatuRn(aSwitchList, reduceToSwitchingGenes = FALSE, alpha = 0.05, dIFcutoff = 0.1, diagplots = FALSE)
} else {
    isoformSwitchTestDEXSeq(aSwitchList, reduceToSwitchingGenes = FALSE, alpha = 0.05, dIFcutoff = 0.1)
}
f <- aSwitchList$isoformFeatures
sig <- f[!is.na(f$isoform_switch_q_value) & f$isoform_switch_q_value < 0.05 & abs(f$dIF) > 0.1, ]
cat(sprintf('\nSignificant switching isoforms (q < 0.05, |dIF| > 0.1): %d in %d genes\n', nrow(sig), length(unique(sig$gene_id))))
if (!demo) cat('Before trusting these calls, run the label-permutation check (SKILL.md "Null check").\n')
write.csv(sig[!duplicated(sig$isoform_id), c('gene_id', 'gene_name', 'isoform_id', 'IF1', 'IF2', 'dIF', 'isoform_switch_q_value')],
          file.path(out_dir, 'significant_switches.csv'), row.names = FALSE)
if (nrow(sig) == 0) {
    cat('No switches at these cutoffs. Nothing to annotate. Compare calculateCountsFromAbundance = TRUE/FALSE (SKILL.md "Count route").\n')
    quit(save = 'no', status = 0)
}

# ---- 3. ORFs: annotated CDS if the GTF has any, else predict --------------------------------------
# analyzeORF() OVERWRITES annotated ORFs (orf_origin -> Predicted) and is less accurate on annotated transcripts.
has_cds <- any(aSwitchList$orfAnalysis$orf_origin == 'Annotation', na.rm = TRUE)
if (!has_cds) aSwitchList <- analyzeORF(aSwitchList, orfMethod = 'longest', genomeObject = NULL, showProgress = FALSE)
cat('ORF origin:', paste(names(table(aSwitchList$orfAnalysis$orf_origin)), table(aSwitchList$orfAnalysis$orf_origin), collapse = ', '), '\n')

# ---- 4. splicing events and sequences for the external annotators --------------------------------
aSwitchList <- analyzeAlternativeSplicing(aSwitchList, onlySwitchingGenes = TRUE, showProgress = FALSE)
seq_dir <- file.path(out_dir, 'sequences'); dir.create(seq_dir, showWarnings = FALSE)
aSwitchList <- extractSequence(aSwitchList, onlySwitchingGenes = TRUE, pathToOutput = seq_dir, writeToFile = TRUE)
cat('Sequences for CPC2 / hmmscan / SignalP / IUPred2A written to', seq_dir, '\n')

# ---- 5. import whichever annotator results exist; request only those consequence types -----------
types <- c('intron_retention', 'ORF_seq_similarity', 'NMD_status')
have <- function(f) !is.null(annot_dir) && file.exists(file.path(annot_dir, f))
if (have('cpc2_result.txt')) {
    # FALSE keeps CPC2-"noncoding" ORFs, which otherwise drops PTC-bearing (NMD) isoforms from the NMD call
    aSwitchList <- analyzeCPC2(aSwitchList, pathToCPC2resultFile = file.path(annot_dir, 'cpc2_result.txt'), removeNoncodinORFs = FALSE)
    types <- c(types, 'coding_potential')
}
if (have('pfam_scanfmt.txt')) {
    aSwitchList <- analyzePFAM(aSwitchList, pathToPFAMresultFile = file.path(annot_dir, 'pfam_scanfmt.txt'))
    types <- c(types, 'domains_identified')
}
if (have('signalp_results.txt')) {
    aSwitchList <- analyzeSignalP(aSwitchList, pathToSignalPresultFile = file.path(annot_dir, 'signalp_results.txt'))
    types <- c(types, 'signal_peptide_identified')
}
if (have('iupred2a_result.txt')) {
    aSwitchList <- analyzeIUPred2A(aSwitchList, pathToIUPred2AresultFile = file.path(annot_dir, 'iupred2a_result.txt'))
    types <- c(types, 'IDR_identified', 'IDR_type')
}
cat('Consequence types requested:', paste(types, collapse = ', '), '\n')
aSwitchList <- analyzeSwitchConsequences(aSwitchList, consequencesToAnalyze = types, dIFcutoff = 0.1, showProgress = FALSE)

# ---- 6. summarise and plot -----------------------------------------------------------------------
print(extractSwitchSummary(aSwitchList, filterForConsequences = TRUE))
top <- extractTopSwitches(aSwitchList, filterForConsequences = FALSE, n = 5, sortByQvals = TRUE)
print(top[, intersect(c('gene_name', 'condition_1', 'condition_2', 'gene_switch_q_value', 'switchConsequencesGene'), colnames(top))])
cond <- levels(factor(design$condition))
plot_file <- file.path(out_dir, sprintf('switchPlot_%s.pdf', top$gene_name[1]))
pdf(plot_file, width = 10, height = 8)
switchPlot(aSwitchList, gene = top$gene_name[1], condition1 = cond[1], condition2 = cond[2])
dev.off()
cat('Plot:', plot_file, '\n')
saveRDS(aSwitchList, file.path(out_dir, 'switchAnalyzeRlist.rds'))

# ---- demo self-check ------------------------------------------------------------------------------
if (demo) {
    found <- intersect(unique(sig$gene_id), toy$planted)
    cat(sprintf('\nDEMO CHECK: planted switching genes recovered %d/%d; other genes called %d\n',
                length(found), length(toy$planted), length(setdiff(unique(sig$gene_id), toy$planted))))
    stopifnot(length(found) >= 12, file.exists(plot_file), file.size(plot_file) > 5000)
    cc <- aSwitchList$switchConsequence
    nmd <- cc[which(grepl('NMD', cc$featureCompared) & cc$switchConsequence == 'NMD sensitive'), ]
    nmd_genes <- unique(sub('_[ABC]$', '', c(nmd$isoformUpregulated, nmd$isoformDownregulated)))
    cat(sprintf('NMD-sensitive switches: %d genes, %d of them planted poison genes (8 planted)\n',
                length(nmd_genes), length(intersect(nmd_genes, toy$poison))))
    stopifnot(setequal(nmd_genes, toy$poison))
    cat('DEMO OK\n')
}
