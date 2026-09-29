#!/usr/bin/env Rscript
# Tested contract: BaalChIP 1.36.0 (Bioconductor 3.22; R 4.5).
# This single-group runner requires either measured per-variant RAF or matching
# gDNA BAM input. It never derives RAF from AF or a detached CN-segment BED.

script_arg <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", script_arg, value = TRUE)
script_path <- if (length(script_flag)) sub("^--file=", "", script_flag[[1]]) else "scripts/baalchip_workflow.R"
source(file.path(dirname(normalizePath(script_path, mustWork = TRUE)), "baalchip_contracts.R"))

usage <- function() {
    cat(paste(
        "Usage: Rscript baalchip_workflow.R",
        "--samples SAMPLE.tsv --hets HETS.tsv --group GROUP",
        "--blacklist BLACKLIST.bed --imprinted IMPRINTED.bed --sex female|male|unknown",
        "--assembly BUILD --samples-assembly BUILD --hets-assembly BUILD",
        "--blacklist-assembly BUILD --imprinted-assembly BUILD",
        "--correction raf|gdna [--gdna-bam GDNA.bam] --samtools /path/to/samtools",
        "--out NEW_OUTPUT_DIR",
        "[--cores N] [--preflight-only]",
        sep = "\n  "
    ), "\n")
}

parse_cli <- function(args) {
    flags <- c("--preflight-only")
    result <- list(preflight_only = FALSE)
    i <- 1L
    while (i <= length(args)) {
        key <- args[[i]]
        if (key %in% flags) {
            result[[gsub("-", "_", sub("^--", "", key))]] <- TRUE
            i <- i + 1L
        } else {
            if (!startsWith(key, "--") || i == length(args)) asb_stop("invalid_argument", paste("invalid CLI token:", key))
            name <- gsub("-", "_", sub("^--", "", key))
            result[[name]] <- args[[i + 1L]]
            i <- i + 2L
        }
    }
    result
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
required <- c(
    "samples", "hets", "group", "blacklist", "imprinted", "sex", "assembly",
    "samples_assembly", "hets_assembly", "blacklist_assembly", "imprinted_assembly",
    "correction", "samtools", "out"
)
missing_args <- required[!required %in% names(args)]
if (length(missing_args)) {
    usage()
    asb_stop("missing_argument", paste("missing:", paste(paste0("--", gsub("_", "-", missing_args)), collapse = ", ")))
}

if (!args$sex %in% c("female", "male", "unknown")) asb_stop("invalid_argument", "--sex must be female, male, or unknown")
if (!args$correction %in% c("raf", "gdna")) asb_stop("invalid_correction", "--correction must be raf or gdna")
if (args$correction == "gdna" && !"gdna_bam" %in% names(args)) asb_stop("missing_gdna", "gDNA correction requires --gdna-bam")
if (args$correction == "raf" && "gdna_bam" %in% names(args)) asb_stop("ambiguous_correction", "choose RAF or gDNA, not both")

validate_assembly(args$assembly, list(
    samples = args$samples_assembly, hets = args$hets_assembly,
    blacklist = args$blacklist_assembly, imprinted = args$imprinted_assembly
))
samples <- validate_samplesheet(args$samples, args$group)
hets <- validate_hets(args$hets, args$correction)
blacklist <- assert_readable_file(args$blacklist, "blacklist BED")
imprinted <- assert_readable_file(args$imprinted, "imprinted-locus BED")

styles <- c(
    hets = hets$contig_style,
    blacklist = contig_style(read_bed_contigs(blacklist, "blacklist BED"), "blacklist BED"),
    imprinted = contig_style(read_bed_contigs(imprinted, "imprinted-locus BED"), "imprinted-locus BED")
)
for (i in seq_along(samples$bed_paths)) {
    styles[[paste0("peaks_", i)]] <- contig_style(
        read_bed_contigs(samples$bed_paths[[i]], paste("peak BED", i)), paste("peak BED", i)
    )
}
for (i in seq_along(samples$bam_paths)) {
    styles[[paste0("bam_", i)]] <- bam_contig_style(samples$bam_paths[[i]], args$samtools, paste("ChIP BAM", i))
}
if (length(unique(styles)) != 1L) {
    asb_stop("contig_mismatch", paste("contig styles disagree:", paste(names(styles), styles, sep = "=", collapse = ", ")))
}

gdna <- list()
if (args$correction == "gdna") {
    gdna_path <- assert_readable_file(args$gdna_bam, "gDNA BAM")
    assert_readable_file(paste0(gdna_path, ".bai"), "gDNA BAM index")
    gdna <- setNames(list(gdna_path), args$group)
    styles[["gdna_bam"]] <- bam_contig_style(gdna_path, args$samtools, "gDNA BAM")
}

if (length(unique(styles)) != 1L) {
    asb_stop("contig_mismatch", paste("contig styles disagree:", paste(names(styles), styles, sep = "=", collapse = ", ")))
}

out <- normalizePath(args$out, winslash = "/", mustWork = FALSE)
if (dir.exists(out) || file.exists(out)) asb_stop("output_exists", paste("refusing to overwrite", out))
parent <- dirname(out)
if (!dir.exists(parent)) asb_stop("missing_output_parent", paste("output parent does not exist:", parent))

run_staged <- function(args, samples, hets, styles, gdna, out, parent) {
stage <- tempfile(pattern = paste0(basename(out), ".partial-"), tmpdir = parent)
if (!dir.create(stage)) asb_stop("write_failure", paste("cannot create staging directory", stage))
published <- FALSE
on.exit(if (!published && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)

status <- data.frame(
    state = "preflight_passed", assembly = args$assembly, group = args$group,
    correction = args$correction, baalchip_version = BAALCHIP_CONTRACT_VERSION, stringsAsFactors = FALSE
)
atomic_write_table(status, file.path(stage, "run-status.tsv"))
atomic_write_table(data.frame(input = names(styles), contig_style = unname(styles)), file.path(stage, "contig-contract.tsv"))
atomic_write_table(data.frame(
    correction = args$correction,
    source = if (args$correction == "raf") hets$path else gdna[[args$group]],
    RAFcorrection = TRUE,
    RMcorrection = TRUE,
    stringsAsFactors = FALSE
), file.path(stage, "correction-provenance.tsv"))

if (isTRUE(args$preflight_only)) {
    if (!file.rename(stage, out)) asb_stop("write_failure", paste("cannot publish", out))
    published <- TRUE
    quit(save = "no", status = 0L)
}

suppressPackageStartupMessages({
    library(BaalChIP)
    library(GenomicRanges)
    library(rtracklayer)
})
assert_baalchip_version(packageVersion("BaalChIP"))

hets_gr <- GRanges(seqnames = hets$table$CHROM, ranges = IRanges(start = hets$table$POS, width = 1))
imprinted_gr <- import(imprinted)
blacklist_gr <- import(blacklist)
reason <- rep(NA_character_, nrow(hets$table))
reason[overlapsAny(hets_gr, imprinted_gr)] <- "imprinted_locus"
reason[is.na(reason) & overlapsAny(hets_gr, blacklist_gr)] <- "blacklist"
if (args$sex == "female") {
    x_contig <- if (identical(hets$contig_style, "chr")) "chrX" else "X"
    reason[is.na(reason) & hets$table$CHROM == x_contig] <- "chrX_female_autosomal_analysis"
}
excluded <- hets$table[!is.na(reason), , drop = FALSE]
excluded$exclusion_reason <- reason[!is.na(reason)]
included <- hets$table[is.na(reason), , drop = FALSE]
atomic_write_table(excluded, file.path(stage, "excluded-variants.tsv"))
filtered_hets <- file.path(stage, "included-variants.tsv")
atomic_write_table(included, filtered_hets)

if (nrow(included) == 0L) {
    status$state <- "no_variants_after_filters"
    atomic_write_table(status, file.path(stage, "run-status.tsv"))
    if (!file.rename(stage, out)) asb_stop("write_failure", paste("cannot publish", out))
    published <- TRUE
    quit(save = "no", status = 3L)
}

res <- BaalChIP(
    samplesheet = samples$path,
    hets = setNames(filtered_hets, args$group),
    CorrectWithgDNA = gdna
)
res <- alleleCounts(res, min_base_quality = 10, min_mapq = 15)
res <- mergePerGroup(res)
res <- filter1allele(res)
merged <- BaalChIP.get(res, "mergedCounts")
if (is.null(merged[[args$group]]) || nrow(merged[[args$group]]) == 0L) {
    empty <- select_group_report(NULL, args$group)
    atomic_write_table(empty, file.path(stage, "ASB-full.tsv"))
    atomic_write_table(empty, file.path(stage, "ASB-significant.tsv"))
    status$state <- "no_calls_after_count_filters"
    status$tested <- 0L
    status$significant <- 0L
    atomic_write_table(status, file.path(stage, "run-status.tsv"))
    if (!file.rename(stage, out)) asb_stop("write_failure", paste("cannot publish", out))
    published <- TRUE
    quit(save = "no", status = 0L)
}
cores <- if ("cores" %in% names(args)) suppressWarnings(as.integer(args$cores)) else 1L
if (is.na(cores) || cores < 1L) asb_stop("invalid_argument", "--cores must be a positive integer")
res <- getASB(
    res, Iter = 5000, conf_level = 0.95, cores = cores,
    RMcorrection = TRUE, RAFcorrection = TRUE
)

asb_table <- select_group_report(BaalChIP.report(res), args$group)
asb_sig <- asb_table[asb_table$isASB %in% TRUE, , drop = FALSE]
atomic_write_table(asb_table, file.path(stage, "ASB-full.tsv"))
atomic_write_table(asb_sig, file.path(stage, "ASB-significant.tsv"))

status$state <- if (nrow(asb_table) == 0L) "no_calls" else "complete"
status$tested <- nrow(asb_table)
status$significant <- nrow(asb_sig)
atomic_write_table(status, file.path(stage, "run-status.tsv"))
if (!file.rename(stage, out)) asb_stop("write_failure", paste("cannot publish", out))
published <- TRUE
cat(sprintf("BaalChIP %s: tested=%d significant=%d correction=%s\n", args$group, nrow(asb_table), nrow(asb_sig), args$correction))
}

run_staged(args, samples, hets, styles, gdna, out, parent)
