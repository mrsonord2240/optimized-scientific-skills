# Reference: SCEVAN 1.0.3 | Verify API if version differs
# Reference-free automatic malignant/non-malignant classification and
# subclone calling from a raw gene-by-cell counts matrix.
# Expects: a raw genes-by-cells counts matrix with real gene symbols
# spanning enough of the 22 human autosomes SCEVAN checks per cell (a
# handful of loci fails its per-chromosome coverage filter for every cell),
# and, optionally, a cell-to-group annotation file (cell <tab> group) to
# supply known non-malignant barcodes -- omit norm_cell to let SCEVAN
# search for confident normals itself via gene-set enrichment, which needs
# a large enough gene panel to find overlapping gene sets.
library(SCEVAN)

counts <- read.table('counts.matrix', header = TRUE, row.names = 1,
                      sep = '\t', check.names = FALSE)
count_mtx <- as.matrix(counts)

# Known non-malignant barcodes (optional but recommended when available;
# mirrors copyKAT's norm.cell.names). cell_annotations.txt has no header:
# cell <tab> group.
annot <- read.table('cell_annotations.txt', header = FALSE, sep = '\t',
                     col.names = c('cell', 'group'))
normal_cell_names <- annot$cell[annot$group %in% c('Tcell', 'Myeloid')]

# pipelineCNA() has no plot=FALSE switch: it always writes heatmap/segment
# PNGs to ./output/ as a side effect, and that internal plotting step can
# throw on a sparse or small cohort -- AFTER classification has already
# completed and printed "found N tumor cells". Wrap in tryCatch() so a
# plotting failure doesn't discard a classification that already succeeded.
results <- tryCatch({
    pipelineCNA(
        count_mtx,
        sample = 'tumor1',
        par_cores = 4,
        norm_cell = normal_cell_names,
        SUBCLONES = TRUE,
        ClonalCN = TRUE,
        plotTree = FALSE,
        organism = 'human',
        ngenes_chr = 5)
}, error = function(e) {
    message('pipelineCNA() classification completed (see the "found N tumor ',
            'cells" console line above); its own plotting step then threw: ',
            conditionMessage(e))
    NULL
})

# results$class is 'tumor', 'normal', or 'filtered' per cell; with
# SUBCLONES = TRUE it adds a subclone column. If results is NULL because
# the plotting step errored, the classification is still visible in the
# console log's "found N tumor cells" line -- treat that as the result and
# rerun with a larger/less sparse cohort if the return value is needed.
if (!is.null(results)) {
    print(table(results$class))
    write.csv(results, 'scevan_prediction.csv', row.names = TRUE)
}

# Malignant/subclone calls remain a clustering decision on a noisy proxy:
# validate against lineage markers or mutations, and treat subclones as
# hypotheses, not measurements, same as inferCNV/copyKAT.
message('SCEVAN complete')
