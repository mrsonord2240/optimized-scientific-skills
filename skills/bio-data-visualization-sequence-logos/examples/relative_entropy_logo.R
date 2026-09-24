# Helpers for a background-corrected DNA/RNA logo in ggseqlogo.
# ggseqlogo's `bits` method uses a uniform background; it does not accept
# a background-frequency argument. These functions compute a custom matrix
# whose column total is D_KL(p || q) in bits and whose letter heights are p_i D_KL.
# It deliberately uses raw observed probabilities: no pseudocount or smoothing
# is applied implicitly.

read_aligned_fasta <- function(path) {
    sequences <- readLines(path, warn = FALSE)
    sequences <- toupper(sequences[nzchar(sequences) & !startsWith(sequences, ">")])
    if (!length(sequences)) stop("No sequences found in ", path, call. = FALSE)
    if (length(unique(nchar(sequences))) != 1L) {
        stop("Sequences must be aligned to one common length", call. = FALSE)
    }
    sequences
}

sequence_probability_matrix <- function(sequences, alphabet) {
    sequences <- toupper(sequences)
    alphabet <- toupper(alphabet)
    if (!length(sequences) || length(unique(nchar(sequences))) != 1L) {
        stop("Sequences must be non-empty and aligned", call. = FALSE)
    }
    observed <- unique(unlist(strsplit(sequences, "", fixed = TRUE)))
    unknown <- setdiff(observed, alphabet)
    if (length(unknown)) {
        stop("Sequences contain symbols absent from the background: ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }

    width <- nchar(sequences[[1]])
    probabilities <- vapply(
        seq_len(width),
        function(position) tabulate(match(substr(sequences, position, position), alphabet),
                                    nbins = length(alphabet)) / length(sequences),
        numeric(length(alphabet))
    )
    rownames(probabilities) <- alphabet
    probabilities
}

relative_entropy_heights <- function(probabilities, background) {
    if (!is.matrix(probabilities) || is.null(rownames(probabilities))) {
        stop("probabilities must be a letter-by-position matrix with row names", call. = FALSE)
    }
    if (is.null(names(background))) {
        stop("background must be a named vector matching the matrix row names", call. = FALSE)
    }
    if (anyDuplicated(rownames(probabilities)) || anyDuplicated(names(background))) {
        stop("probability and background letter names must be unique", call. = FALSE)
    }
    if (!setequal(rownames(probabilities), names(background))) {
        stop("background letter names must exactly match probability row names", call. = FALSE)
    }
    background <- background[rownames(probabilities)]
    if (anyNA(background) || any(!is.finite(background)) || any(background <= 0)) {
        stop("background must give every letter a finite positive frequency", call. = FALSE)
    }
    if (any(!is.finite(probabilities)) || any(probabilities < 0) ||
        any(colSums(probabilities) <= 0)) {
        stop("probabilities must contain non-negative, non-empty columns", call. = FALSE)
    }

    probabilities <- sweep(probabilities, 2, colSums(probabilities), "/")
    background <- background / sum(background)
    contributions <- matrix(0, nrow(probabilities), ncol(probabilities),
                            dimnames = dimnames(probabilities))
    present <- probabilities > 0
    contributions[present] <- probabilities[present] *
        log2(probabilities[present] / background[row(contributions)[present]])
    relative_entropy <- colSums(contributions)
    heights <- sweep(probabilities, 2, relative_entropy, "*")
    attr(heights, "relative_entropy_bits") <- relative_entropy
    heights
}
