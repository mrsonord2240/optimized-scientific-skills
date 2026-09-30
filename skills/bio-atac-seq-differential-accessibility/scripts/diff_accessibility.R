#!/usr/bin/env Rscript
# Reference: DiffBind 3.12+, DESeq2 1.42+, edgeR 4.0+, ChIPseeker 1.38+, sva 3.50+ | Verify API if version differs
# DiffBind workflow with summits=250 fixed-width counting and a >= min_overlap consensus.
# Usage: Rscript diff_accessibility.R <sample_sheet.csv> [output_prefix] [normalize] [fdr_thr] [lfc_thr] [--key=value ...]
#   normalize: lib (default, full library size) | native | tmm | rle | default (DiffBind's own; RiP library)
#   --library=full|peaks    library-size definition for normalize=lib (default full)
#   --method=deseq2|edger   DiffBind backend (default deseq2)
#   --case=treated --control=control   Condition labels in the sheet (default treated / control)
#   --design='~Tissue + Condition'     DiffBind design; terms limited to Tissue, Factor, Condition, Treatment, Replicate
#                                      sheet columns (put batch/donor in one of them)
#   --sva=N                 up to N surrogate variables (capped by sample count); fits DESeq2 directly on the
#                           DiffBind counts as ~SV1+..+Condition (ignores --method/--design; no blacklist filter)
#   --txdb=<TxDb package>|none   annotation source (default TxDb.Hsapiens.UCSC.hg38.knownGene)
#   --min-overlap=2 --summit=250 --cores=1
# The sample sheet needs SampleID, Condition, Replicate, bamReads, Peaks, PeakCaller columns.

suppressPackageStartupMessages({ library(DiffBind); library(rtracklayer); library(ChIPseeker) })

NORM <- list(lib = DBA_NORM_LIB, native = DBA_NORM_NATIVE, tmm = DBA_NORM_TMM,
             rle = DBA_NORM_RLE, default = DBA_NORM_DEFAULT)

run_diff <- function(sample_sheet, output_prefix = 'diff_atac',
                     normalize = 'lib', library_size = 'full',
                     fdr_thr = 0.05, lfc_thr = 1,           # ENCODE-style thresholds
                     method = 'deseq2', case = 'treated', control = 'control',
                     design = NULL, n_sv = 0,
                     min_overlap = 2,                       # peak in >= 2 replicates
                     fixed_width_summit = 250,              # +/- 250 = 501 bp (Corces 2018)
                     txdb = 'TxDb.Hsapiens.UCSC.hg38.knownGene', cores = 1) {
    stopifnot(normalize %in% names(NORM), library_size %in% c('full', 'peaks'),
              method %in% c('deseq2', 'edger'))
    be <- if (method == 'edger') DBA_EDGER else DBA_DESEQ2
    cat('=== DiffBind workflow ===\n')
    dba <- dba(sampleSheet = sample_sheet)
    meta <- dba.show(dba)
    if (!all(c(control, case) %in% meta$Condition))
        stop(sprintf("Condition labels '%s'/'%s' not found; sheet has: %s (set --case/--control)",
                     case, control, paste(unique(meta$Condition), collapse = ', ')))

    cat(sprintf('Counting reads in fixed-width consensus (summits +/- %d, minOverlap=%d)...\n',
                fixed_width_summit, min_overlap))
    dba <- dba.count(dba, summits = fixed_width_summit, minOverlap = min_overlap, bParallel = cores > 1)

    # lib = library-size scaling (full library keeps global shifts); native = RLE (DESeq2) / TMM (edgeR).
    lib_arg <- if (normalize == 'lib') (if (library_size == 'peaks') DBA_LIBSIZE_PEAKREADS else DBA_LIBSIZE_FULL)
               else DBA_LIBSIZE_DEFAULT
    dba <- dba.normalize(dba, normalize = NORM[[normalize]], library = lib_arg)
    nrm <- dba.normalize(dba, bRetrieve = TRUE)
    cat(sprintf('Normalization: norm.method=%s lib.method=%s\n', nrm$norm.method, nrm$lib.method))

    if (n_sv > 0) {
        # DiffBind designs only accept its categorical metadata columns, so surrogate variables
        # go through a direct DESeq2 fit on the DiffBind consensus counts.
        res_sva <- fit_sva(dba, nrm, meta, n_sv, design, case, control, fdr_thr, lfc_thr)
        results <- res_sva$results
    } else {
        cat('Setting up contrast...\n')
        dba <- dba.contrast(dba, design = if (is.null(design)) TRUE else design,
                            contrast = c('Condition', case, control), minMembers = 2)
        cat(sprintf('Fitting %s model...\n', method))
        dba <- dba.analyze(dba, method = be)
        cat(sprintf('Reporting at FDR < %.2f, |log2FC| >= %.1f\n', fdr_thr, lfc_thr))
        results <- dba.report(dba, method = be, th = fdr_thr, fold = lfc_thr)
    }
    n_hit <- if (is.null(results)) 0L else length(results)
    cat(sprintf('  Differentially accessible peaks: %d\n', n_hit))

    pdf(sprintf('%s_diagnostics.pdf', output_prefix), width = 8, height = 6)
    dba.plotPCA(dba, attributes = DBA_CONDITION, label = DBA_ID)
    if (n_sv > 0) {
        plot_ma_volcano(res_sva$table, fdr_thr, lfc_thr)
    } else {
        dba.plotMA(dba, method = be, th = fdr_thr, fold = lfc_thr)
        dba.plotVolcano(dba, method = be, th = fdr_thr, fold = lfc_thr)
        if (n_hit >= 2) dba.plotHeatmap(dba, contrast = 1, method = be, th = fdr_thr, correlations = FALSE)  # heatmap has no fold filter
    }
    dev.off()

    if (n_hit == 0) {
        warning('No sites pass the thresholds; BED files are empty and no annotated CSV is written. Loosen fdr_thr/lfc_thr or check the design.')
        file.create(sprintf('%s_opened.bed', output_prefix), sprintf('%s_closed.bed', output_prefix))
        return(invisible(list(dba = dba, results = results, anno = NULL)))
    }
    cat(sprintf('  Opened (Fold > 0): %d   Closed (Fold < 0): %d\n',
                sum(results$Fold > 0), sum(results$Fold < 0)))
    export(results[results$Fold > 0], sprintf('%s_opened.bed', output_prefix))
    export(results[results$Fold < 0], sprintf('%s_closed.bed', output_prefix))

    peakAnno <- NULL
    if (!identical(txdb, 'none')) {
        cat('Annotating peaks to genes...\n')
        suppressPackageStartupMessages(library(txdb, character.only = TRUE))
        peakAnno <- annotatePeak(results, TxDb = get(txdb), tssRegion = c(-2000, 500),
                                 level = 'gene', verbose = FALSE)
        write.csv(as.data.frame(peakAnno), sprintf('%s_annotated.csv', output_prefix), row.names = FALSE)
        pdf(sprintf('%s_annoplot.pdf', output_prefix))
        plotAnnoPie(peakAnno); print(plotDistToTSS(peakAnno))
        dev.off()
    } else {
        write.csv(as.data.frame(results), sprintf('%s_annotated.csv', output_prefix), row.names = FALSE)
    }
    invisible(list(dba = dba, results = results, anno = peakAnno))
}

fit_sva <- function(dba, nrm, meta, n_sv, design, case, control, fdr_thr, lfc_thr) {
    suppressPackageStartupMessages({ library(DESeq2); library(sva) })
    cat('Estimating surrogate variables with svaseq...\n')
    se <- dba(dba, bSummarizedExperiment = TRUE)
    cnt <- round(SummarizedExperiment::assay(se, 'Reads')); colnames(cnt) <- meta$ID
    mod  <- model.matrix(~Condition, meta)
    mod0 <- model.matrix(~1, meta)
    k <- min(n_sv, nrow(meta) - ncol(mod) - 1)             # svaseq returns NaN above n - ncol(mod) - 1
    if (k < 1) stop(sprintf('%d samples leave no room for surrogate variables', nrow(meta)))
    sv <- sva::svaseq(cnt[rowMeans(cnt) > 1, ], mod, mod0, n.sv = k)$sv
    if (anyNA(sv)) stop('svaseq returned NaN surrogate variables; lower --sva')
    coldata <- data.frame(Condition = factor(meta$Condition, levels = c(control, case)), row.names = meta$ID)
    coldata[paste0('SV', seq_len(k))] <- sv
    dsn <- as.formula(paste('~', paste0('SV', seq_len(k), collapse = ' + '), '+ Condition'))
    if (!is.null(design)) message('--design is ignored with --sva; the model is ', deparse(dsn))
    cat(sprintf('  -> %d SV(s); DESeq2 design %s\n', k, deparse(dsn)))
    dds <- DESeqDataSetFromMatrix(cnt, coldata, dsn)
    sizeFactors(dds) <- if (nrm$norm.method == 'lib') nrm$lib.sizes / min(nrm$lib.sizes) else nrm$norm.factors
    dds <- DESeq(dds, quiet = TRUE)
    tab <- as.data.frame(results(dds, contrast = c('Condition', case, control)))
    gr <- SummarizedExperiment::rowRanges(se)
    S4Vectors::mcols(gr) <- data.frame(Conc = log2(tab$baseMean + 1), Fold = tab$log2FoldChange,
                                       p.value = tab$pvalue, FDR = tab$padj)
    keep <- !is.na(tab$padj) & tab$padj < fdr_thr & abs(tab$log2FoldChange) >= lfc_thr
    list(results = if (any(keep)) gr[keep] else NULL, table = tab)
}

plot_ma_volcano <- function(tab, fdr_thr, lfc_thr) {
    sig <- !is.na(tab$padj) & tab$padj < fdr_thr & abs(tab$log2FoldChange) >= lfc_thr
    col <- ifelse(sig, 'red3', 'grey60'); ttl <- sprintf('%d sites FDR<%.2f, |log2FC|>=%.1f', sum(sig), fdr_thr, lfc_thr)
    plot(log2(tab$baseMean + 1), tab$log2FoldChange, pch = 20, col = col, xlab = 'log2 mean count',
         ylab = 'log2FC', main = paste('MA:', ttl)); abline(h = 0)
    plot(tab$log2FoldChange, -log10(tab$pvalue), pch = 20, col = col, xlab = 'log2FC',
         ylab = '-log10 p', main = paste('Volcano:', ttl))
}

parse_cli <- function(args) {
    flags <- grep('^--', args, value = TRUE); pos <- grep('^--', args, value = TRUE, invert = TRUE)
    opt <- setNames(sub('^--[^=]+=', '', flags), sub('^--([^=]+)=.*$', '\\1', flags))
    unknown <- setdiff(names(opt), c('library', 'method', 'case', 'control', 'design', 'sva',
                                     'txdb', 'min-overlap', 'summit', 'cores'))
    if (length(unknown) || length(pos) < 1 || length(pos) > 5)
        stop('Usage: diff_accessibility.R <sample_sheet.csv> [prefix] [normalize] [fdr] [lfc] [--key=value ...]; bad: ',
             paste(c(unknown, if (length(pos) > 5) 'extra positional'), collapse = ', '))
    get <- function(k, d) if (k %in% names(opt)) opt[[k]] else d
    a <- list(sample_sheet = pos[1])
    if (length(pos) >= 2) a$output_prefix <- pos[2]
    if (length(pos) >= 3) a$normalize <- tolower(sub('^DBA_NORM_', '', toupper(pos[3])))
    if (length(pos) >= 4) a$fdr_thr <- as.numeric(pos[4])
    if (length(pos) >= 5) a$lfc_thr <- as.numeric(pos[5])
    a$library_size <- get('library', 'full'); a$method <- get('method', 'deseq2')
    a$case <- get('case', 'treated'); a$control <- get('control', 'control')
    if ('design' %in% names(opt)) a$design <- opt[['design']]
    a$n_sv <- as.integer(get('sva', 0)); a$txdb <- get('txdb', 'TxDb.Hsapiens.UCSC.hg38.knownGene')
    a$min_overlap <- as.integer(get('min-overlap', 2)); a$fixed_width_summit <- as.integer(get('summit', 250))
    a$cores <- as.integer(get('cores', 1))
    if (anyNA(c(a$fdr_thr, a$lfc_thr))) stop('fdr_thr and lfc_thr must be numeric')
    a
}

if (sys.nframe() == 0) {                       # run only as a script; source() to reuse run_diff()
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) do.call(run_diff, parse_cli(args))
}
