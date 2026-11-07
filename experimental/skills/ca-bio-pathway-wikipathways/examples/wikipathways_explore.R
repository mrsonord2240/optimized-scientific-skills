# Reference: clusterProfiler 4.18+, rWikiPathways 1.26+ | Verify API if version differs
# NOTE: every rWikiPathways query and download hits the live WikiPathways API/archive over the
# network (needs internet). downloadPathwayArchive pins a dated MONTHLY release = the reproducible
# pattern; enrichWP/gseWP without a pin pull current/, which changes month to month.
library(rWikiPathways)
library(clusterProfiler)
library(tidyr)

organisms <- rWikiPathways::listOrganisms()  # supported species; choose its exact scientific-name string
head(organisms)

human_pathways <- listPathways('Homo sapiens')
head(human_pathways, 20)

# searchPathways() is NOT a current function (errors); use findPathwaysByText
cancer_pathways <- findPathwaysByText('cancer')
head(cancer_pathways)

pathway_info <- getPathwayInfo('WP554')   # metadata incl. last-edit; check before trusting a hit
pathway_info

# genes by BridgeDb system code: 'L'=Entrez, 'H'=HGNC symbol, 'En'=Ensembl
pathway_entrez <- getXrefList('WP554', 'L')
pathway_entrez

# Reproducible pattern: pin a dated GMT, split the compound term, run enricher on the pinned sets.
# format defaults to gpml -> pass format='gmt'; organism=NULL would open the index, not download.
# Releases land on the 10th of each month and the live archive keeps only ~12 months (see
# data.wikipathways.org's own landing page) -- compute a recent date instead of hardcoding one
# that will 404 later; for a fixed historical date beyond the window, use the Zenodo GMT/GPML
# archive instead (https://zenodo.org/communities/wikipathways).
# Try newest-first in-window dates; one missing release 404s, so step back a month and retry.
candidates <- unique(format(Sys.Date() - seq(60, 330, by = 30), '%Y%m10'))
gmt <- NULL
for (archive_date in candidates) {   # report the date that succeeds in methods
  gmt <- tryCatch(suppressWarnings(downloadPathwayArchive(date = archive_date, organism = 'Homo sapiens',
                                                          format = 'gmt', destpath = tempdir())),
                  error = function(e) NULL)
  if (!is.null(gmt) && file.exists(file.path(tempdir(), gmt))) break
  gmt <- NULL
}
if (is.null(gmt)) stop('no release in the last ~12 months; use the Zenodo GMT archive')
wp2gene <- read.gmt(file.path(tempdir(), gmt))
wp2gene <- separate(wp2gene, term, c('name', 'version', 'wpid', 'org'), sep = '%')
t2g <- wp2gene[, c('wpid', 'gene')]
t2n <- wp2gene[, c('wpid', 'name')]

# self-contained demo of the pinned-archive pattern: reuse the WP554 genes fetched above (line 23)
# as the query set against the pinned archive's own gene universe (a real analysis instead passes
# the tested-gene universe -- see wikipathways_ora.R); report date=archive_date in methods
entrez_ids <- pathway_entrez
all_entrez <- unique(t2g$gene)
wp_pinned <- enricher(entrez_ids, universe = all_entrez, TERM2GENE = t2g, TERM2NAME = t2n)
as.data.frame(wp_pinned)
