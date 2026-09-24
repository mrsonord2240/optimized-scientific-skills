# Purpose: mummichog / PSEA on a raw m/z peak table (no compound IDs) with MetaboAnalystR. The peak
#          table MUST be the ENTIRE feature table (R_all), never significant features only.
# Inputs:  peak table (columns m/z, p-value, t-score; 'mpt' format); optional ppm (5.0), ionization
#          mode (negative|positive, default negative), query p-cutoff (0.2), permNum (1000),
#          output CSV (mummichog_psea.csv). Needs fitdistrplus and RJSONIO installed.
# Output:  the CSV of predicted-active pathways (mSet$mummi.resmat); NOT a metabolite ID list.
# Usage:   Rscript scripts/mummichog_psea.R peaks.csv [5.0] [negative] [0.2] [1000] [mummichog_psea.csv]
# Note:    downloads the mummichog library (e.g. hsa_mfn.qs) into the working directory on first run.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop('usage: Rscript mummichog_psea.R peaks.csv [ppm] [mode] [pval] [permNum] [out]')
peaks_file <- args[1]
ppm      <- if (length(args) >= 2) as.numeric(args[2]) else 5.0
mode     <- if (length(args) >= 3) args[3] else 'negative'
pcut     <- if (length(args) >= 4) as.numeric(args[4]) else 0.2
perm_num <- if (length(args) >= 5) as.integer(args[5]) else 1000L
out_file <- if (length(args) >= 6) args[6] else 'mummichog_psea.csv'

current.msg <- character(0); err.vec <- character(0)  # required -- see Version Compatibility
library(MetaboAnalystR)
set.seed(123)  # PerformPSEA's permNum resampling is not reproducible run-to-run otherwise

mSet <- InitDataObjects('mass_all', 'mummichog', FALSE)
mSet <- SetPeakFormat(mSet, 'mpt')               # 'mpt' = m/z, p-value, t-score; 'mprt' adds RT (use with 'v2')

# ppm and ionization mode are chemistry-specific and mandatory; pos and neg use
# entirely different adduct tables. Mixed data needs a per-feature mode column.
mSet <- UpdateInstrumentParameters(mSet, ppm, mode)

# CRITICAL: peaks_file must be the ENTIRE feature table, not just significant peaks.
# The permutation null draws random feature lists from this file (R_all); supplying
# only significant features pre-enriches the pool and makes everything significant.
mSet <- Read.PeakListData(mSet, peaks_file)
mSet <- SanityCheckMummichogData(mSet)          # also merges duplicate m/z-matched features automatically

mSet <- SetPeakEnrichMethod(mSet, 'mum', 'v2')   # 'mum'|'gsea'|'integ'; 'v2' uses RT/empirical compounds
mSet <- SetMummichogPval(mSet, pcut)             # query-defining cutoff; default is NOT 0.05 -- document it
mSet <- PerformPSEA(mSet, 'hsa_mfn', 'current', permNum = perm_num) # library string encodes organism+network

psea <- mSet$mummi.resmat                         # predicted-active pathways; NOT a metabolite ID list
write.csv(as.data.frame(psea), out_file)
print(head(as.data.frame(psea), 10))
