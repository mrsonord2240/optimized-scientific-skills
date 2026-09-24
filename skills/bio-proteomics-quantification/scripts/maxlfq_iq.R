# Run the real MaxLFQ (iq) on a long peptide-ion table; per-run median normalization first.
# Input  : CSV with columns protein, ion, run, intensity (RAW intensities, one row per peptide ion per run)
# Output : CSV protein x sample matrix of log2 MaxLFQ estimates, plus <out>.disconnected.txt listing proteins
#          whose samples are NOT on one common scale (non-empty create_protein_table annotation)
# Usage  : Rscript scripts/maxlfq_iq.R peptide_long.csv protein_maxlfq.csv
# Checked: iq 2.0.1 (R 4.4.3)
suppressPackageStartupMessages(library(iq))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('usage: maxlfq_iq.R peptide_long.csv out.csv')
peptide_long <- read.csv(args[1], stringsAsFactors = FALSE)

# peptide_long: one row per peptide ion per run, RAW intensities (protein, ion, run, intensity)
norm <- preprocess(peptide_long, primary_id = 'protein', secondary_id = 'ion', sample_id = 'run',
                   intensity_col = 'intensity', median_normalization = TRUE, pdf_out = NULL)
protein_list  <- create_protein_list(norm)    # one run-normalized log2 matrix per protein
protein_table <- create_protein_table(protein_list, method = 'maxLFQ')
protein_matrix <- protein_table$estimate      # proteins x samples

# maxLFQ solves per-protein least squares only within CONNECTED sample sets; a non-empty annotation
# marks proteins whose samples split into groups that are NOT on one common scale
disconnected <- rownames(protein_matrix)[nzchar(protein_table$annotation)]

write.csv(protein_matrix, args[2])
writeLines(disconnected, paste0(args[2], '.disconnected.txt'))
cat('proteins:', nrow(protein_matrix), '| samples:', ncol(protein_matrix), '| disconnected:', length(disconnected), '\n')
