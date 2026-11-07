#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 1L) stop("usage: test_baalchip_contracts.R PATH_TO_CONTRACTS", call. = FALSE)
source(argv[[1]])

expect_error_code <- function(expr, code) {
    found <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
    if (is.null(found) || !grepl(sprintf("^[[]%s[]]", code), found)) {
        stop(sprintf("expected [%s], got: %s", code, if (is.null(found)) "success" else found), call. = FALSE)
    }
}

stopifnot(
    identical(BAALCHIP_CONTRACT_VERSION, "1.36.0"),
    identical(BAALCHIP_CONTRACT_BIOCONDUCTOR, "3.22"),
    isTRUE(assert_baalchip_version("1.36.0"))
)
expect_error_code(assert_baalchip_version("1.38.0"), "version_mismatch")

root <- tempfile("baal-contracts-")
dir.create(root)
on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

writeLines("dummy", file.path(root, "r1.bam"))
writeLines("dummy", file.path(root, "r1.bam.bai"))
writeLines("chr1\t0\t100", file.path(root, "peaks.bed"))
samples <- data.frame(
    group_name = "TUMOR", target = "FOXA1", replicate_number = 1L,
    bam_name = "r1.bam", bed_name = "peaks.bed", stringsAsFactors = FALSE
)
write.table(samples, file.path(root, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
valid_samples <- validate_samplesheet(file.path(root, "samples.tsv"), "TUMOR")
stopifnot(nrow(valid_samples$table) == 1L)
expect_error_code(validate_samplesheet(file.path(root, "samples.tsv"), "OTHER"), "sample_mismatch")

hets <- data.frame(ID = "v1", CHROM = "chr1", POS = 10L, REF = "A", ALT = "G", RAF = 0.7)
write.table(hets, file.path(root, "hets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
valid_raf <- validate_hets(file.path(root, "hets.tsv"), "raf")
stopifnot(valid_raf$contig_style == "chr", valid_raf$table$RAF == 0.7)
valid_gdna <- validate_hets(file.path(root, "hets.tsv"), "gdna")
stopifnot(!"RAF" %in% colnames(valid_gdna$table))

names(hets)[names(hets) == "RAF"] <- "AF"
write.table(hets, file.path(root, "af-not-raf.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
expect_error_code(validate_hets(file.path(root, "af-not-raf.tsv"), "raf"), "missing_raf")
expect_error_code(validate_assembly("GRCh38", list(samples = "GRCh38", hets = "hg19")), "assembly_mismatch")

report_columns <- c(
    "ID", "CHROM", "POS", "REF", "ALT", "REF.counts", "ALT.counts",
    "Total.counts", "AR", "RMbias", "RAF", "Bayes_lower", "Bayes_upper",
    "Corrected.AR", "isASB"
)
row <- as.data.frame(setNames(as.list(seq_along(report_columns)), report_columns), check.names = FALSE)
selected <- select_group_report(list(TUMOR = row), "TUMOR")
stopifnot(nrow(selected) == 1L, identical(colnames(selected), report_columns))
stopifnot(nrow(select_group_report(NULL, "TUMOR")) == 0L)
expect_error_code(select_group_report(row, "TUMOR"), "report_contract")

cat("baalchip contract tests: PASS\n")
