#!/usr/bin/env Rscript
# Purpose: per-cell motif deviations and per-cluster marker motifs in an ArchR project.
# Usage:   Rscript chromvar_archr.R project_dir markers.rds [group=Clusters] [motifset=cisbp] [macs2=<path>] [peakset=existing|build]
#            peakset=build runs addReproduciblePeakSet(groupBy=group, pathToMacs2=macs2) first (needs macs2=).
#            Default peakset=existing: the project already has a peak set.
# Input:   saved ArchR project directory (loadArchRProject) with cell groups in the `group` column.
# Output:  markers.rds = getMarkerFeatures(useMatrix='MotifMatrix', useSeqnames='z'); the project gains
#          PeakMatrix, Motif annotations, background peaks and MotifMatrix. Cells with NA z-scores are
#          dropped from the project in memory only and counted on stdout.
# Checked: ArchR 1.0.3 (not yet executed after the recut).
suppressPackageStartupMessages(library(ArchR))

args <- commandArgs(trailingOnly = TRUE)
opts <- list(group = 'Clusters', motifset = 'cisbp', macs2 = '', peakset = 'existing')
positional <- character(0)
for (a in args) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) >= 2 && kv[1] %in% names(opts)) opts[[kv[1]]] <- paste(kv[-1], collapse = '=')
  else if (grepl('=', a, fixed = TRUE)) stop('unknown argument ', a)
  else positional <- c(positional, a)
}
if (length(positional) != 2) stop('usage: chromvar_archr.R project_dir markers.rds [group=Clusters] [motifset=cisbp] [macs2=path] [peakset=existing|build]')
proj <- loadArchRProject(positional[1])
if (!opts$group %in% names(getCellColData(proj))) stop('project has no cell column "', opts$group, '"')

if (opts$peakset == 'build') {
  if (!nzchar(opts$macs2)) stop('peakset=build needs macs2=<path to macs2>')
  proj <- addReproduciblePeakSet(proj, groupBy = opts$group, pathToMacs2 = opts$macs2)
}
proj <- addPeakMatrix(proj)  # use the PeakMatrix; tile deviations add intergenic noise
proj <- addMotifAnnotations(proj, motifSet = opts$motifset, name = 'Motif')
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(proj, peakAnnotation = 'Motif')

# Sparse cells x rare motifs can give 0/0 (NaN) z-scores, and getMarkerFeatures refuses NA.
mm <- getMatrixFromProject(proj, useMatrix = 'MotifMatrix')
na_cells <- colnames(mm)[colSums(is.na(assays(mm)$z)) > 0]
cat('cells with NA z-scores:', length(na_cells), 'of', ncol(mm), '\n')
if (length(na_cells) > 0) proj <- subsetCells(proj, cellNames = setdiff(proj$cellNames, na_cells))

markers <- getMarkerFeatures(proj, useMatrix = 'MotifMatrix', groupBy = opts$group, useSeqnames = 'z')
saveRDS(markers, positional[2])
message('wrote ', positional[2])
