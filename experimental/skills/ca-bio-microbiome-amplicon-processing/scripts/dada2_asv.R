#!/usr/bin/env Rscript
# Checked on DADA2 1.30+
# Per-run DADA2 ASV inference: filterAndTrim -> learnErrors -> dada -> mergePairs, then
# mergeSequenceTables across runs and ONE removeBimeraDenovo.
# Primers MUST already be removed (scripts/remove_primers.sh; ITSxpress for ITS).
#
# Usage: Rscript dada2_asv.R --trunc-len F,R [options] run1_dir [run2_dir ...]
#   --trunc-len F,R   required for paired reads (truncLen is a merge budget, see routes/truncation.md).
#                     With --single: one number; 0 for ITS (never fixed-truncate ITS)
#   --single          single-end input (ITSxpress output, one FASTQ per sample); no pair merging
#   --max-ee 2,2      maxEE (a single number with --single); loosen the reverse for long amplicons
#   --trunc-q 2       truncQ
#   --min-len 0       minLen (50 for ITS)
#   --amplicon-len N  if given, refuse a paired truncLen that cannot overlap (sum < N + 12)
#   --pool FALSE      dada() pool: FALSE, TRUE or pseudo (pseudo for rare-ASV sensitivity)
#   --pattern-f/--pattern-r  FASTQ suffixes (defaults _R1_001.fastq.gz / _R2_001.fastq.gz; single default .fastq.gz)
#   --out dada2_out   output directory
# Each run_dir holds the primer-trimmed FASTQs of ONE sequencing run; the error model is fit per run.
# Writes seqtab_nochim.rds, read_tracking.csv, error_fit_*.png, asv_lengths.txt in --out.

suppressPackageStartupMessages({
  library(dada2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list(`trunc-len` = NULL, single = FALSE, `max-ee` = NULL, `trunc-q` = "2", `min-len` = "0",
            `amplicon-len` = NULL, pool = "FALSE", `pattern-f` = NULL, `pattern-r` = "_R2_001.fastq.gz",
            out = "dada2_out")
runs <- character()
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a == "--single") {
    opt$single <- TRUE
  } else if (startsWith(a, "--")) {
    key <- sub("^--", "", a)
    if (!(key %in% names(opt)) || i == length(args)) stop("unknown or incomplete option: ", a)
    opt[[key]] <- args[i + 1]
    i <- i + 1
  } else {
    runs <- c(runs, a)
  }
  i <- i + 1
}
if (is.null(opt$`pattern-f`)) opt$`pattern-f` <- if (opt$single) ".fastq.gz" else "_R1_001.fastq.gz"
if (length(runs) == 0 || is.null(opt$`trunc-len`)) {
  stop("Usage: Rscript dada2_asv.R --trunc-len F,R [--single] [--max-ee 2,2] [--amplicon-len N] [--pool pseudo] run1_dir [run2_dir ...]")
}

num <- function(x) as.numeric(strsplit(x, ",")[[1]])
truncLen <- num(opt$`trunc-len`)
maxEE <- if (is.null(opt$`max-ee`)) rep(2, if (opt$single) 1 else 2) else num(opt$`max-ee`)
truncQ <- as.numeric(opt$`trunc-q`)
minLen <- as.numeric(opt$`min-len`)
pool <- switch(opt$pool, `FALSE` = FALSE, `TRUE` = TRUE, pseudo = "pseudo", stop("--pool must be FALSE, TRUE or pseudo"))

if (!opt$single && length(truncLen) != 2) stop("paired input needs --trunc-len F,R")
if (opt$single && length(truncLen) != 1) stop("--single needs one --trunc-len value (0 for ITS)")
if (!opt$single && !is.null(opt$`amplicon-len`)) {
  need <- as.numeric(opt$`amplicon-len`) + 12
  if (sum(truncLen) < need) {
    stop(sprintf("truncLen %d + %d = %d is below amplicon length + 12 = %d: pairs cannot merge. Keep more length, loosen maxEE on the reverse.",
                 truncLen[1], truncLen[2], sum(truncLen), need))
  }
}
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

sample_names_of <- function(fn) {
  sapply(strsplit(sub("\\.fastq(\\.gz)?$", "", basename(fn)), "_"), `[`, 1)
}

# A read shorter than truncLen is dropped, not padded: truncLen above the post-primer read length empties the filter.
check_ceiling <- function(fn, trunc, label) {
  if (all(trunc == 0)) return(invisible())
  fs <- ShortRead::FastqSampler(fn[1], 1000)
  w <- nchar(as.character(ShortRead::sread(ShortRead::yield(fs))))
  close(fs)
  message(sprintf("%s: primer-trimmed read length in %s: min %d, median %d, max %d", label, basename(fn[1]), min(w), median(w), max(w)))
  if (trunc > max(w)) stop(sprintf("%s truncLen %d exceeds the longest primer-trimmed read (%d): every read would be dropped", label, trunc, max(w)))
}

process_run <- function(run_path) {
  fnFs <- sort(list.files(run_path, pattern = paste0(opt$`pattern-f`, "$"), full.names = TRUE))
  if (length(fnFs) == 0) stop("no files matching ", opt$`pattern-f`, " in ", run_path)
  sample_names <- sample_names_of(sub(opt$`pattern-f`, ".fastq.gz", fnFs, fixed = TRUE))
  filtFs <- file.path(run_path, "filtered", paste0(sample_names, "_F.fastq.gz"))
  names(filtFs) <- sample_names
  run_name <- basename(normalizePath(run_path))

  if (opt$single) {
    check_ceiling(fnFs, truncLen[1], "reads")
    out <- filterAndTrim(fnFs, filtFs, truncLen = truncLen, maxEE = maxEE, truncQ = truncQ,
                         minLen = minLen, maxN = 0, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
    errF <- learnErrors(filtFs, multithread = TRUE)
    ggsave(file.path(opt$out, paste0("error_fit_", run_name, ".png")), plotErrors(errF, nominalQ = TRUE), width = 8, height = 6)
    dadaFs <- dada(filtFs, err = errF, pool = pool, multithread = TRUE)
    seqtab <- makeSequenceTable(dadaFs)
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN))
    colnames(track) <- c("input", "filtered", "denoised")
  } else {
    fnRs <- sort(list.files(run_path, pattern = paste0(opt$`pattern-r`, "$"), full.names = TRUE))
    if (length(fnRs) != length(fnFs)) stop("forward and reverse file counts differ in ", run_path)
    filtRs <- file.path(run_path, "filtered", paste0(sample_names, "_R.fastq.gz"))
    names(filtRs) <- sample_names
    check_ceiling(fnFs, truncLen[1], "forward")
    check_ceiling(fnRs, truncLen[2], "reverse")
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = truncLen, maxEE = maxEE, truncQ = truncQ,
                         minLen = minLen, maxN = 0, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
    errF <- learnErrors(filtFs, multithread = TRUE)   # fit THIS run only
    errR <- learnErrors(filtRs, multithread = TRUE)
    # Observed points must track the fitted line and fall with Q; on binned-quality NovaSeq/NextSeq data they may not.
    ggsave(file.path(opt$out, paste0("error_fit_F_", run_name, ".png")), plotErrors(errF, nominalQ = TRUE), width = 8, height = 6)
    ggsave(file.path(opt$out, paste0("error_fit_R_", run_name, ".png")), plotErrors(errR, nominalQ = TRUE), width = 8, height = 6)
    dadaFs <- dada(filtFs, err = errF, pool = pool, multithread = TRUE)
    dadaRs <- dada(filtRs, err = errR, pool = pool, multithread = TRUE)
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
    seqtab <- makeSequenceTable(mergers)
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
  }
  rownames(track) <- sample_names
  list(seqtab = seqtab, track = track)
}

per_run <- lapply(runs, process_run)

# Combine AFTER per-run inference; the exact sequence string is the join key.
seqtab_all <- if (length(per_run) > 1) mergeSequenceTables(tables = lapply(per_run, `[[`, "seqtab")) else per_run[[1]]$seqtab

# ONE chimera removal on the combined table. Chimeras are many ASVs but few reads; a large read loss is a leftover-primer smell.
seqtab_nochim <- removeBimeraDenovo(seqtab_all, method = "consensus", multithread = TRUE, verbose = TRUE)
cat("reads retained after chimera removal:", round(100 * sum(seqtab_nochim) / sum(seqtab_all), 1), "%\n")

lens <- table(nchar(getSequences(seqtab_nochim)))   # off-target lengths flag mis-merges or non-target amplification
capture.output(print(lens), file = file.path(opt$out, "asv_lengths.txt"))
print(lens)

write.csv(do.call(rbind, lapply(per_run, `[[`, "track")), file.path(opt$out, "read_tracking.csv"))
saveRDS(seqtab_nochim, file.path(opt$out, "seqtab_nochim.rds"))   # carry 'run' forward as a batch covariate
cat("Done:", file.path(opt$out, "seqtab_nochim.rds"), "\n")
