# Purpose: build a cohort-complete ComplexHeatmap oncoPrint matrix from MAF-like rows.
# Inputs: a data.frame (or tab-separated file), an explicit cohort sample vector, and optional normalized CNV/fusion calls.
# Usage: source("scripts/maf_to_oncoprint.R"); mat <- maf_to_oncoprint(maf, clinical$Tumor_Sample_Barcode)

ONCOPRINT_CLASS_MAP <- c(
  Missense_Mutation = "Missense",
  In_Frame_Ins = "Missense",
  In_Frame_Del = "Missense",
  Nonsense_Mutation = "Truncating",
  Frame_Shift_Ins = "Truncating",
  Frame_Shift_Del = "Truncating",
  Nonstop_Mutation = "Truncating",
  Translation_Start_Site = "Truncating",
  Splice_Site = "Splice"
)

ONCOPRINT_CLASSES <- c("Amp", "HomDel", "Missense", "Truncating", "Splice", "Fusion")

.require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(label, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

maf_to_oncoprint <- function(maf,
                             cohort_samples,
                             genes = NULL,
                             top = 20L,
                             additional_calls = NULL,
                             class_map = ONCOPRINT_CLASS_MAP) {
  if (is.character(maf) && length(maf) == 1L) {
    maf <- read.delim(maf, comment.char = "#", stringsAsFactors = FALSE,
                      check.names = FALSE, quote = "")
  }
  .require_columns(maf,
                   c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"),
                   "maf")

  cohort_samples <- as.character(cohort_samples)
  if (!length(cohort_samples) || anyNA(cohort_samples) || any(cohort_samples == "")) {
    stop("cohort_samples must contain non-empty sample IDs.", call. = FALSE)
  }
  if (anyDuplicated(cohort_samples)) {
    stop("cohort_samples contains duplicate sample IDs.", call. = FALSE)
  }

  mapped <- unname(class_map[as.character(maf$Variant_Classification)])
  ignored <- sort(unique(as.character(maf$Variant_Classification[is.na(mapped)])))
  keep <- !is.na(mapped) & maf$Tumor_Sample_Barcode %in% cohort_samples
  calls <- data.frame(
    Hugo_Symbol = as.character(maf$Hugo_Symbol[keep]),
    Tumor_Sample_Barcode = as.character(maf$Tumor_Sample_Barcode[keep]),
    alteration_class = mapped[keep],
    stringsAsFactors = FALSE
  )

  if (!is.null(additional_calls)) {
    .require_columns(additional_calls,
                     c("Hugo_Symbol", "Tumor_Sample_Barcode", "alteration_class"),
                     "additional_calls")
    additional_calls <- additional_calls[, c(
      "Hugo_Symbol", "Tumor_Sample_Barcode", "alteration_class"
    )]
    additional_calls[] <- lapply(additional_calls, as.character)
    unknown <- setdiff(unique(additional_calls$alteration_class), ONCOPRINT_CLASSES)
    if (length(unknown)) {
      stop("additional_calls contains unsupported alteration class(es): ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    calls <- rbind(calls, additional_calls[
      additional_calls$Tumor_Sample_Barcode %in% cohort_samples, , drop = FALSE
    ])
  }
  calls <- unique(calls)

  if (is.null(genes)) {
    if (!nrow(calls)) stop("No mapped alterations remain for the cohort.", call. = FALSE)
    gene_frequency <- tapply(calls$Tumor_Sample_Barcode, calls$Hugo_Symbol,
                             function(x) length(unique(x)))
    genes <- names(sort(gene_frequency, decreasing = TRUE))
    genes <- head(genes, as.integer(top))
  }
  genes <- unique(as.character(genes))
  if (!length(genes) || anyNA(genes) || any(genes == "")) {
    stop("genes must contain at least one non-empty gene symbol.", call. = FALSE)
  }

  mat <- matrix("", nrow = length(genes), ncol = length(cohort_samples),
                dimnames = list(genes, cohort_samples))
  calls <- calls[calls$Hugo_Symbol %in% genes, , drop = FALSE]
  if (nrow(calls)) {
    grouped <- split(calls$alteration_class,
                     paste(calls$Hugo_Symbol, calls$Tumor_Sample_Barcode, sep = "\r"))
    for (key in names(grouped)) {
      parts <- strsplit(key, "\r", fixed = TRUE)[[1]]
      values <- unique(grouped[[key]])
      values <- ONCOPRINT_CLASSES[ONCOPRINT_CLASSES %in% values]
      mat[parts[1], parts[2]] <- paste(values, collapse = ";")
    }
  }
  if (!any(mat != "")) {
    stop("The selected genes have no mapped alterations; oncoPrint cannot render an all-empty matrix.",
         call. = FALSE)
  }
  attr(mat, "ignored_variant_classifications") <- ignored
  mat
}

align_oncoprint_clinical <- function(clinical,
                                     sample_ids,
                                     sample_col = "Tumor_Sample_Barcode") {
  .require_columns(clinical, sample_col, "clinical")
  ids <- as.character(clinical[[sample_col]])
  if (anyDuplicated(ids)) stop("clinical contains duplicate sample IDs.", call. = FALSE)
  matched <- match(sample_ids, ids)
  if (anyNA(matched)) {
    stop("Clinical metadata is missing matrix sample(s): ",
         paste(head(sample_ids[is.na(matched)], 10L), collapse = ", "), call. = FALSE)
  }
  aligned <- clinical[matched, , drop = FALSE]
  stopifnot(identical(as.character(aligned[[sample_col]]), as.character(sample_ids)))
  rownames(aligned) <- NULL
  aligned
}
