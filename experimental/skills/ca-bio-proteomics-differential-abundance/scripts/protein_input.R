# Shared input helpers for limma_de.R, deqms_de.R and proda_de.R (sourced, not run).
# Reads either a raw search-engine table (MaxQuant proteinGroups-style TSV: linear `Intensity <sample>`
# columns, empty or 0 = missing, `Reverse` / `Potential contaminant` / `Only identified by site` flags)
# or a log2 matrix CSV (first column = protein ID, one column per sample, NA = missing).

MEDIAN_LINEAR <- 50  # a log2 intensity matrix has a median near 15-35; above this the values are linear

parse_kv <- function(args) {
  kv <- list()
  for (a in args) {
    i <- regexpr('=', a, fixed = TRUE)
    if (i < 2) stop('Unrecognised argument: ', a, ' (options are key=value)')
    kv[[substr(a, 1, i - 1)]] <- substr(a, i + 1, nchar(a))
  }
  kv
}

read_samples <- function(path) {
  s <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c('sample', 'condition') %in% names(s))) stop('samples.csv needs columns: sample, condition [, batch]')
  s
}

MIN_COMPLETE <- 50  # complete-case proteins needed before the sample median is taken on them alone

# Subtract each sample's median (location only), taken over proteins observed in EVERY sample when there
# are enough of them: under intensity-dependent missingness a sample that lost more low-abundance proteins
# has an upward-biased median over its own observed set.
center_medians <- function(protein_matrix) {
  complete <- rowSums(is.na(protein_matrix)) == 0
  reference <- if (sum(complete) >= MIN_COMPLETE) protein_matrix[complete, , drop = FALSE] else protein_matrix
  medians <- apply(reference, 2, median, na.rm = TRUE)
  sweep(protein_matrix, 2, medians - mean(medians))
}

# Returns list(matrix = log2 proteins x samples (columns in samples order), counts = named numeric or NULL,
#              removed = data.frame(protein, reason) for flagged rows).
read_protein_matrix <- function(path, samples, kv = list()) {
  format <- if (!is.null(kv$format)) kv$format else if (grepl('\\.(tsv|txt|tab)$', path, ignore.case = TRUE)) 'raw' else 'log2'
  removed <- data.frame(protein = character(0), reason = character(0))
  counts <- NULL
  if (format == 'raw') {
    prefix <- if (!is.null(kv$prefix)) kv$prefix else 'Intensity '
    tab <- read.delim(path, check.names = FALSE, quote = '', comment.char = '', stringsAsFactors = FALSE)
    cols <- paste0(prefix, samples$sample)
    miss <- setdiff(cols, names(tab))
    if (length(miss)) stop('Columns not in the table: ', paste(head(miss, 5), collapse = ', '),
                           '. Intensity-like columns: ',
                           paste(head(grep('ntensity', names(tab), value = TRUE), 5), collapse = ', '),
                           '. Set prefix= (e.g. prefix="LFQ intensity ") so prefix + sample name matches a column.')
    id_col <- if (!is.null(kv$id)) kv$id else intersect(c('Protein IDs', 'Majority protein IDs'), names(tab))[1]
    if (is.na(id_col)) id_col <- names(tab)[1]
    ids <- make.unique(as.character(tab[[id_col]]))
    drop <- rep(FALSE, nrow(tab))
    for (flag in c('Reverse', 'Potential contaminant', 'Only identified by site')) {
      if (flag %in% names(tab)) {
        hit <- tab[[flag]] %in% '+' & !drop  # empty cells are NA; '+' marks the row
        removed <- rbind(removed, data.frame(protein = ids[hit], reason = rep(flag, sum(hit))))
        drop <- drop | hit
      }
    }
    m <- sapply(cols, function(cl) suppressWarnings(as.numeric(tab[[cl]])))
    m <- matrix(m, nrow = nrow(tab), dimnames = list(ids, samples$sample))
    m[!is.na(m) & m <= 0] <- NA  # 0 is undetected, not a measurement
    m <- log2(m)
    if (!is.null(kv$count_col)) {
      if (!kv$count_col %in% names(tab)) stop('count_col not in the table: ', kv$count_col)
      counts <- setNames(suppressWarnings(as.numeric(tab[[kv$count_col]])), ids)[!drop]
    }
    m <- m[!drop, , drop = FALSE]
  } else {
    m <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
    miss <- setdiff(samples$sample, colnames(m))
    if (length(miss)) stop('Samples not in the matrix columns: ', paste(head(miss, 5), collapse = ', '))
    m <- m[, samples$sample, drop = FALSE]
    mid <- median(m, na.rm = TRUE)
    if (mid > MEDIAN_LINEAR) stop(sprintf(paste(
      'Median value is %.0f: these look like linear intensities, not log2. Give the raw search-engine TSV',
      '(.tsv with Intensity columns), or take log2 first.'), mid))
  }
  list(matrix = m, counts = counts, removed = removed)
}
