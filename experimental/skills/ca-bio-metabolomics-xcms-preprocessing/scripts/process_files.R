args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('usage: process_files.R sample_metadata.csv output-directory')
library(xcms)
library(MsExperiment)
pd <- read.csv(args[1], stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c('file_path', 'sample_name', 'sample_group') %in% names(pd)) ||
    anyNA(pd[, c('file_path', 'sample_name', 'sample_group')]) ||
    anyDuplicated(pd$sample_name) || !all(file.exists(pd$file_path))) stop('invalid sample metadata or paths')
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
raw <- readMsExperiment(spectraFiles = pd$file_path, sampleData = pd)
cwp <- CentWaveParam(ppm = 25, peakwidth = c(20, 80), snthresh = 10,
                     prefilter = c(3, 1000), noise = 1000)
xdata <- findChromPeaks(raw, param = cwp)
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
pdp <- PeakDensityParam(sampleGroups = sampleData(xdata)$sample_group,
                        bw = 30, minFraction = 0.5, binSize = 0.025)
xdata <- groupChromPeaks(xdata, param = pdp)
write.csv(featureValues(xdata, value = 'into'), file.path(args[2], 'unfilled.csv'))
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())
write.csv(featureValues(xdata, value = 'into'), file.path(args[2], 'filled.csv'))
write.csv(as.data.frame(featureDefinitions(xdata)), file.path(args[2], 'features.csv'))
saveRDS(xdata, file.path(args[2], 'xcms.rds'))
