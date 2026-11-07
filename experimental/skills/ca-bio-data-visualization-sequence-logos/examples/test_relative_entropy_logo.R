# Deterministic numerical and ggseqlogo API checks for relative_entropy_logo.R.
# Run from this examples directory.

source("relative_entropy_logo.R")

background <- c(A = 0.29, C = 0.21, G = 0.21, T = 0.29)
probabilities <- cbind(
    pure_c = c(A = 0, C = 1, G = 0, T = 0),
    background_matched = background
)
heights <- relative_entropy_heights(probabilities, background)

stopifnot(
    isTRUE(all.equal(sum(heights[, "pure_c"]), -log2(background[["C"]]), tolerance = 1e-12)),
    isTRUE(all.equal(sum(heights[, "background_matched"]), 0, tolerance = 1e-12)),
    isTRUE(all.equal(attr(heights, "relative_entropy_bits"), colSums(heights), tolerance = 1e-12)),
    all(is.finite(heights)),
    all(heights[c("A", "G", "T"), "pure_c"] == 0)
)

mixed <- cbind(mixed = c(A = 0.50, C = 0.25, G = 0.25, T = 0))
mixed_heights <- relative_entropy_heights(mixed, background)
mixed_relative_entropy <- sum(mixed[, "mixed"][mixed[, "mixed"] > 0] *
                              log2(mixed[, "mixed"][mixed[, "mixed"] > 0] /
                                   background[mixed[, "mixed"] > 0]))
stopifnot(isTRUE(all.equal(
    mixed_heights[mixed[, "mixed"] > 0, "mixed"],
    mixed[mixed[, "mixed"] > 0, "mixed"] * mixed_relative_entropy,
    tolerance = 1e-12
)))

uniform_background <- setNames(rep(0.25, 4), names(background))
uniform_heights <- relative_entropy_heights(mixed, uniform_background)
uniform_relative_entropy <- 2 + sum(mixed[, "mixed"][mixed[, "mixed"] > 0] *
                                    log2(mixed[, "mixed"][mixed[, "mixed"] > 0]))
stopifnot(isTRUE(all.equal(
    uniform_heights[mixed[, "mixed"] > 0, "mixed"],
    mixed[mixed[, "mixed"] > 0, "mixed"] * uniform_relative_entropy,
    tolerance = 1e-12
)))

gc_rich_background <- c(A = 0.18, C = 0.32, G = 0.32, T = 0.18)
gc_calibration <- relative_entropy_heights(cbind(
    pure_g = c(A = 0, C = 0, G = 1, T = 0),
    pure_a = c(A = 1, C = 0, G = 0, T = 0)
), gc_rich_background)
stopifnot(
    isTRUE(all.equal(sum(gc_calibration[, "pure_g"]), 1.643856189774725, tolerance = 1e-12)),
    isTRUE(all.equal(sum(gc_calibration[, "pure_a"]), 2.473931188332412, tolerance = 1e-12)),
    isTRUE(all.equal(
        relative_entropy_heights(mixed, background * 10),
        mixed_heights, tolerance = 1e-12
    )),
    inherits(try(relative_entropy_heights(mixed, c(A = .5, C = .5)), silent = TRUE), "try-error"),
    inherits(try(relative_entropy_heights(mixed, c(A = .2, A = .2, C = .2, G = .2)), silent = TRUE), "try-error"),
    inherits(try(relative_entropy_heights(mixed, c(background, N = .01)), silent = TRUE), "try-error"),
    inherits(try(relative_entropy_heights(rbind(mixed, N = 0), background), silent = TRUE), "try-error")
)

fixture <- read_aligned_fasta("aligned_motif.fa")
fixture_heights <- relative_entropy_heights(
    sequence_probability_matrix(fixture, names(background)), background
)
stopifnot(
    nrow(fixture_heights) == 4L,
    ncol(fixture_heights) == nchar(fixture[[1]]),
    identical(rownames(fixture_heights), names(background)),
    all(is.finite(fixture_heights)),
    all(colSums(fixture_heights) >= 0)
)

library(ggseqlogo)
custom_plot <- ggseqlogo(heights[, "pure_c", drop = FALSE], method = "custom")
plotted_layer <- ggplot2::ggplot_build(custom_plot)$data[[1]]
stopifnot(
    inherits(custom_plot, "ggplot"),
    isTRUE(all.equal(max(plotted_layer$y), -log2(background[["C"]]), tolerance = 1e-12))
)
message("relative-entropy logo checks passed")
