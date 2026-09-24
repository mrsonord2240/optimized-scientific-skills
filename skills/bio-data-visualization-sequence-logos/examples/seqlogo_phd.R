# Reference: ggseqlogo 0.2+ | Run from this examples directory.
#
# This example distinguishes conventional uniform-background bits from a
# background-corrected relative-entropy logo. ggseqlogo has no background
# frequency argument: for a non-uniform background, compute letter heights
# explicitly and pass them with method = "custom".

library(ggseqlogo)
library(ggplot2)
library(patchwork)
source("relative_entropy_logo.R")

human_bg <- c(A = 0.29, C = 0.21, G = 0.21, T = 0.29)
dna_alphabet <- names(human_bg)

custom_relative_entropy_logo <- function(sequences) {
    relative_entropy_heights(
        sequence_probability_matrix(sequences, dna_alphabet),
        human_bg
    )
}

# 1. INPUT -- bundled, aligned DNA sequences
seqs <- read_aligned_fasta("aligned_motif.fa")
n <- length(seqs)

# 2. BACKGROUND-CORRECTED LOGO
# The column total is D_KL(observed composition || human_bg), in bits.
p_relative_entropy <- ggseqlogo(custom_relative_entropy_logo(seqs),
                                method = "custom",
                                col_scheme = "nucleotide") +
    labs(title = sprintf("Motif logo (relative entropy; N = %d; human background)", n),
         y = "bits relative to background") +
    theme(plot.title = element_text(size = 10))

# 3. UNIFORM-BACKGROUND BITS LOGO FOR COMPARISON
p_bits <- ggseqlogo(seqs, method = "bits") +
    labs(title = sprintf("Same motif (uniform-background bits; N = %d)", n),
         y = "bits") +
    theme(plot.title = element_text(size = 10))

ggsave("logo_comparison.pdf", p_relative_entropy / p_bits,
       width = 89, height = 80, units = "mm", device = cairo_pdf)

# 4. MULTI-MOTIF STACK -- labels are derived from the bundled fixture data.
motif_sequences <- list(
    CTCF = read_aligned_fasta("ctcf_aligned.fa"),
    REST = read_aligned_fasta("rest_aligned.fa"),
    GATA1 = read_aligned_fasta("gata1_aligned.fa")
)
multi <- lapply(motif_sequences, custom_relative_entropy_logo)
names(multi) <- sprintf("%s (N=%d)", names(motif_sequences),
                        vapply(motif_sequences, length, integer(1)))
p_stack <- ggseqlogo(multi, method = "custom", col_scheme = "nucleotide", ncol = 1) +
    labs(y = "bits relative to human background")
ggsave("multi_logo.pdf", p_stack, width = 89, height = 110, units = "mm",
       device = cairo_pdf)

# 5. PROTEIN LOGO with a custom functional-class palette.
phospho_neighborhoods <- read_aligned_fasta("phospho_aligned.fa")
protein_scheme <- make_col_scheme(
    chars = c("S", "T", "Y",                         # phospho-acceptors
              "K", "R", "H",                         # basic
              "D", "E",                              # acidic
              "A", "V", "L", "I", "M", "F", "W", "C", "G", "P", "N", "Q"),
    cols = c(rep("#D55E00", 3), rep("#0072B2", 3), rep("#CC79A7", 2),
             rep("#009E73", 12)))
p_protein <- ggseqlogo(phospho_neighborhoods, method = "bits", seq_type = "aa",
                        col_scheme = protein_scheme) +
    labs(title = sprintf("Phospho-substrate neighborhoods (N = %d)",
                         length(phospho_neighborhoods)))
ggsave("protein_logo.pdf", p_protein, width = 130, height = 50, units = "mm",
       device = cairo_pdf)

# 6. PWM INPUT: rows are letters, columns are positions. Convert counts to
# probabilities before calling relative_entropy_heights for a corrected logo.
pwm_counts <- matrix(c(85, 5, 5, 5, 5, 85, 5, 5, 30, 20, 30, 20, 5, 5, 5, 85),
                     ncol = 4, byrow = FALSE,
                     dimnames = list(c("A", "C", "G", "T"), NULL))
pwm_probabilities <- sweep(pwm_counts, 2, colSums(pwm_counts), "/")
p_pwm <- ggseqlogo(relative_entropy_heights(pwm_probabilities, human_bg), method = "custom",
                   col_scheme = "nucleotide") +
    labs(y = "bits relative to human background")
ggsave("pwm_logo.pdf", p_pwm, width = 89, height = 50, units = "mm",
       device = cairo_pdf)
