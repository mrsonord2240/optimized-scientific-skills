# Stage 1 -- xcms 4.x feature extraction: centroided mzML -> feat (features x samples) + defs.
# Inputs (define before source()): mzml_files (character vector of centroided mzML paths),
#   pd (data.frame, one row per file in the same order, with a sample_group column 'QC'/'Control'/...),
#   ionization_mode ('positive' or 'negative').
# Outputs (left in the calling environment): xdata, feat, defs.
# Usage: source('scripts/stage1_xcms_extract.R')
# Checked on xcms 4.4.0 (R 4.4.3).
stopifnot(exists('mzml_files'), exists('pd'), exists('ionization_mode'),
          length(mzml_files) == nrow(pd), 'sample_group' %in% colnames(pd))
library(xcms); library(MsExperiment)   # xcms 4.4.0 does not attach MsExperiment: readMsExperiment()/sampleData() are missing without it
# pd: data.frame, one row per file, with a sample_group column ('QC'/'Control'/'Treatment')
# ionization_mode: 'positive' or 'negative', set from the acquisition method -- carried through
# to Stage 3/5 as defs$mode below (commitment #1, the mode-lock). A mixed-mode study runs this
# whole stage twice, once per mode, and merges the resulting feature tables afterward.
raw <- readMsExperiment(spectraFiles = mzml_files, sampleData = pd)

cwp <- CentWaveParam(ppm = 10, peakwidth = c(2, 20), snthresh = 10,
                     prefilter = c(3, 1000), noise = 1000)   # set from instrument; see xcms-preprocessing
xdata <- findChromPeaks(raw, param = cwp)
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6,
    subset = which(sampleData(xdata)$sample_group == 'QC'), subsetAdjust = 'average'))   # anchor RT alignment on pooled QCs
pdp <- PeakDensityParam(sampleGroups = sampleData(xdata)$sample_group,
                        bw = 5, minFraction = 0.5, binSize = 0.025)
xdata <- groupChromPeaks(xdata, param = pdp)        # group on corrected RT (obiwarp needs no pre-grouping)
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())

feat <- featureValues(xdata, value = 'into')        # features x samples; filled cells are imputations
defs <- featureDefinitions(xdata)                   # mzmed / rtmed per feature, for annotation + mummichog
