# tmt_impurity_correct.R -- TMT10 reporter extraction from mzML and lot-specific CoA impurity correction.
# Purpose : the single-plex TMT route of bio-workflows-proteomics-pipeline, with the CoA orientation guard.
# Inputs  : mzML with MS2 reporter ions; tmt10_coa.csv (square percentage matrix, rows = SOURCE reagent,
#           columns = OBSERVED channel, both named 126, 127N, 127C, ...)
# Usage   : Rscript tmt_impurity_correct.R [tmt.mzML] [tmt10_coa.csv] [reporters.csv]
# Output  : reporters.csv = purity-corrected reporter intensities per spectrum
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
mzml_file <- arg(1, 'tmt.mzML')
coa_file  <- arg(2, 'tmt10_coa.csv')
out_file  <- arg(3, 'tmt_reporters_corrected.csv')

library(MSnbase)

# Extract reporter ions from spectra (NOT readMSnSet, which loads an existing text matrix)
raw <- readMSData(mzml_file, mode = 'onDisk')
tmt_data <- quantify(raw, reporters = TMT10, method = 'max')
# Correct isobaric impurity cross-talk with the LOT-SPECIFIC matrix from the reagent CoA.
# makeImpuritiesMatrix(x = 10, edit = FALSE) returns a MANUFACTURER TEMPLATE, not an identity
# matrix -- its diagonal runs 0.928-0.965 and 5% of 126 lands in 127C. It is a shape check only;
# the numbers are lot-specific. (edit = TRUE, the default, opens an editor and blocks in scripts.)
# Do NOT use makeImpuritiesMatrix(filename = ...) for TMT10/TMTpro: it reads a CoA laid out by
# Da OFFSET and places each column k POSITIONS away in the reporter list, which is only correct
# for non-interleaved reagents (TMT6, iTRAQ). TMT10/TMTpro interleave N and C, so the +1 Da
# neighbour of 126 is 127C -- TWO positions away -- and the filename route silently writes the
# bleed into 127N instead. Build the matrix by CHANNEL NAME and hand it to purityCorrect:
# coa_file: a square percentage matrix, rows = SOURCE reagent, columns = OBSERVED channel,
# both labelled with the channel names (126, 127N, 127C, ...); diagonal = the lot's purity.
coa <- as.matrix(read.csv(coa_file, row.names = 1, check.names = FALSE))
stopifnot(nrow(coa) == ncol(coa), setequal(rownames(coa), reporterNames(TMT10)))
coa <- coa[reporterNames(TMT10), reporterNames(TMT10)]   # force the quant's channel order
# ORIENTATION CHECK. A transposed sheet has the right channel names, the right shape and produces
# NO negative values, so the negative-count check below never sees it -- yet it makes the
# correction 2.7x worse than the correct orientation (median relative error 0.0015 vs 0.0006 on
# the TMT10 template), still better than doing nothing and therefore silent. Test the orientation
# directly: a ROW is one reagent's isotopic envelope and sums to 100% by construction (minus what
# falls off the ends of the channel list); a COLUMN sums over different reagents and has no such
# constraint. If the columns fit 100 better than the rows, the sheet is the wrong way round.
row_dev <- sum((rowSums(coa) - 100)^2); col_dev <- sum((colSums(coa) - 100)^2)
if (col_dev < row_dev)
    stop('CoA looks TRANSPOSED: column sums fit 100% better than row sums (', round(col_dev, 1),
         ' vs ', round(row_dev, 1), '). Rows must be the SOURCE reagent, columns the OBSERVED ',
         'channel -- transpose the sheet or re-read the lot certificate.')
stopifnot(all(diag(coa) == apply(coa, 1, max)),   # each reagent's own channel must dominate its row
          all(diag(coa) > 50))                    # a CoA is percentages; < 50 means fractions were read
impurities <- coa / 100                                  # CoA percentages -> fractions
tmt_data <- purityCorrect(tmt_data, impurities)
stopifnot(sum(exprs(tmt_data) < 0, na.rm = TRUE) == 0)   # negatives = a grossly wrong matrix (NOT a transposition test; see above)
write.csv(exprs(tmt_data), out_file)
cat('wrote', out_file, ':', nrow(tmt_data), 'spectra x', ncol(tmt_data), 'channels\n')
