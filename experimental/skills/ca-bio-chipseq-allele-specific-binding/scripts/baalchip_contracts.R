asb_stop <- function(code, message) {
    stop(sprintf("[%s] %s", code, message), call. = FALSE)
}

BAALCHIP_CONTRACT_VERSION <- "1.36.0"
BAALCHIP_CONTRACT_BIOCONDUCTOR <- "3.22"

assert_baalchip_version <- function(installed_version) {
    installed <- as.character(installed_version)
    if (!identical(installed, BAALCHIP_CONTRACT_VERSION)) {
        asb_stop("version_mismatch", sprintf(
            "this runner is contract-tested against BaalChIP %s (Bioconductor %s); found %s",
            BAALCHIP_CONTRACT_VERSION, BAALCHIP_CONTRACT_BIOCONDUCTOR, installed
        ))
    }
    invisible(TRUE)
}

assert_scalar <- function(value, name) {
    if (length(value) != 1L || is.na(value) || !nzchar(value)) {
        asb_stop("invalid_argument", sprintf("%s must be one non-empty value", name))
    }
    invisible(value)
}

assert_readable_file <- function(path, name) {
    assert_scalar(path, name)
    if (!file.exists(path) || file.info(path)$isdir) {
        asb_stop("missing_input", sprintf("%s does not exist: %s", name, path))
    }
    normalizePath(path, winslash = "/", mustWork = TRUE)
}

resolve_from <- function(path, directory) {
    if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(directory, path)
}

assert_columns <- function(tab, required, label) {
    missing <- setdiff(required, colnames(tab))
    if (length(missing) > 0L) {
        asb_stop("schema_mismatch", sprintf(
            "%s is missing required columns: %s", label, paste(missing, collapse = ",")
        ))
    }
    invisible(tab)
}

contig_style <- function(contigs, label) {
    values <- unique(as.character(contigs[!is.na(contigs) & nzchar(contigs)]))
    if (length(values) == 0L) asb_stop("empty_contigs", sprintf("%s has no contigs", label))
    has_chr <- grepl("^chr", values)
    if (any(has_chr) && !all(has_chr)) {
        asb_stop("mixed_contigs", sprintf("%s mixes chr-prefixed and unprefixed contigs", label))
    }
    if (all(has_chr)) "chr" else "bare"
}

read_bed_contigs <- function(path, label) {
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(lines) & !grepl("^(#|track[[:space:]]|browser[[:space:]])", lines)]
    if (length(lines) == 0L) asb_stop("empty_intervals", sprintf("%s has no data rows", label))
    fields <- strsplit(lines, "[[:space:]]+")
    if (any(lengths(fields) < 3L)) {
        asb_stop("schema_mismatch", sprintf("%s must have at least three BED columns", label))
    }
    vapply(fields, `[[`, character(1), 1L)
}

bam_contig_style <- function(path, samtools, label) {
    samtools <- assert_readable_file(samtools, "samtools executable")
    header <- system2(samtools, c("view", "-H", path), stdout = TRUE, stderr = TRUE)
    status <- attr(header, "status")
    if (!is.null(status) && status != 0L) {
        asb_stop("invalid_bam", sprintf("cannot read %s header", label))
    }
    sq <- header[startsWith(header, "@SQ")]
    names <- sub(".*(?:^|\\t)SN:([^\\t]+).*", "\\1", sq, perl = TRUE)
    if (length(sq) == 0L || any(names == sq)) {
        asb_stop("invalid_bam", sprintf("%s header has no parseable @SQ/SN contigs", label))
    }
    contig_style(names, label)
}

validate_assembly <- function(expected, declarations) {
    assert_scalar(expected, "analysis assembly")
    bad <- names(declarations)[vapply(declarations, function(x) {
        length(x) != 1L || is.na(x) || x != expected
    }, logical(1))]
    if (length(bad) > 0L) {
        asb_stop("assembly_mismatch", sprintf(
            "declared assembly must be %s for: %s", expected, paste(bad, collapse = ",")
        ))
    }
    invisible(expected)
}

validate_samplesheet <- function(path, group) {
    path <- assert_readable_file(path, "samplesheet")
    tab <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("group_name", "target", "replicate_number", "bam_name", "bed_name")
    assert_columns(tab, required, "samplesheet")
    if (nrow(tab) == 0L) asb_stop("empty_samples", "samplesheet has no rows")
    if (any(!nzchar(trimws(as.matrix(tab[, required]))))) {
        asb_stop("schema_mismatch", "samplesheet required fields cannot be empty")
    }
    if (!all(tab$group_name == group)) {
        asb_stop("sample_mismatch", "this single-group workflow requires every group_name to equal --group")
    }
    reps <- suppressWarnings(as.integer(tab$replicate_number))
    if (any(is.na(reps) | reps < 1L) || anyDuplicated(reps)) {
        asb_stop("sample_mismatch", "replicate_number must contain unique positive integers")
    }
    root <- dirname(path)
    bam_paths <- vapply(tab$bam_name, resolve_from, character(1), directory = root)
    bed_paths <- vapply(tab$bed_name, resolve_from, character(1), directory = root)
    for (p in bam_paths) {
        assert_readable_file(p, "ChIP BAM")
        assert_readable_file(paste0(p, ".bai"), "ChIP BAM index")
    }
    for (p in bed_paths) assert_readable_file(p, "peak BED")
    list(path = path, table = tab, bam_paths = bam_paths, bed_paths = bed_paths)
}

validate_hets <- function(path, correction) {
    path <- assert_readable_file(path, "heterozygous variant table")
    tab <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("ID", "CHROM", "POS", "REF", "ALT")
    assert_columns(tab, required, "heterozygous variant table")
    if (nrow(tab) == 0L) asb_stop("empty_variants", "heterozygous variant table has no rows")
    if (anyDuplicated(tab$ID) || any(!nzchar(tab$ID))) {
        asb_stop("schema_mismatch", "variant ID values must be non-empty and unique")
    }
    pos <- suppressWarnings(as.integer(tab$POS))
    if (any(is.na(pos) | pos < 1L)) asb_stop("schema_mismatch", "POS must be positive 1-based integers")
    if (any(!toupper(tab$REF) %in% c("A", "C", "G", "T")) ||
        any(!toupper(tab$ALT) %in% c("A", "C", "G", "T")) ||
        any(nchar(tab$REF) != 1L | nchar(tab$ALT) != 1L)) {
        asb_stop("schema_mismatch", "REF and ALT must be single uppercase A/C/G/T bases")
    }
    if (any(toupper(tab$REF) == toupper(tab$ALT))) {
        asb_stop("schema_mismatch", "REF and ALT must differ")
    }
    if (identical(correction, "raf")) {
        if (!"RAF" %in% colnames(tab)) {
            if ("AF" %in% colnames(tab)) {
                asb_stop("missing_raf", "AF is not RAF; supply measured per-variant reference allele frequency")
            }
            asb_stop("missing_raf", "RAF correction requires a RAF column")
        }
        raf <- suppressWarnings(as.numeric(tab$RAF))
        if (any(!is.finite(raf) | raf < 0 | raf > 1)) {
            asb_stop("invalid_raf", "RAF values must be finite and within [0,1]")
        }
        tab$RAF <- raf
    } else if (identical(correction, "gdna")) {
        # BaalChIP gives a supplied RAF column priority over CorrectWithgDNA.
        # Remove it so the requested gDNA correction cannot be silently bypassed.
        tab$RAF <- NULL
    } else {
        asb_stop("invalid_correction", "--correction must be raf or gdna")
    }
    tab$CHROM <- as.character(tab$CHROM)
    tab$POS <- pos
    tab$REF <- toupper(tab$REF)
    tab$ALT <- toupper(tab$ALT)
    list(path = path, table = tab, contig_style = contig_style(tab$CHROM, "heterozygous variants"))
}

select_group_report <- function(report, group) {
    required <- c(
        "ID", "CHROM", "POS", "REF", "ALT", "REF.counts", "ALT.counts",
        "Total.counts", "AR", "RMbias", "RAF", "Bayes_lower", "Bayes_upper",
        "Corrected.AR", "isASB"
    )
    if (is.null(report)) {
        return(as.data.frame(setNames(replicate(length(required), logical(0), simplify = FALSE), required)))
    }
    if (!is.list(report) || is.data.frame(report) || is.null(names(report)) || !group %in% names(report)) {
        asb_stop("report_contract", "BaalChIP.report must return a named list containing --group")
    }
    result <- report[[group]]
    if (is.null(result)) {
        return(as.data.frame(setNames(replicate(length(required), logical(0), simplify = FALSE), required)))
    }
    if (!is.data.frame(result)) asb_stop("report_contract", "selected group report is not a data.frame")
    assert_columns(result, required, "BaalChIP group report")
    result
}

atomic_write_table <- function(tab, path) {
    tmp <- tempfile(pattern = ".write-", tmpdir = dirname(path))
    on.exit(unlink(tmp), add = TRUE)
    write.table(tab, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    if (!file.rename(tmp, path)) asb_stop("write_failure", sprintf("cannot publish %s", path))
    invisible(path)
}
