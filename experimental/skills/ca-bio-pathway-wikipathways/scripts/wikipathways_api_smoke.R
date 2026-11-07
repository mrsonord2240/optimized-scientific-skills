# Focused live API regression: rWikiPathways standalone organism discovery plus a non-human ORA.
suppressMessages({library(rWikiPathways); library(clusterProfiler)})

organisms <- rWikiPathways::listOrganisms()
stopifnot('Danio rerio' %in% organisms)

zfish_paths <- listPathways('Danio rerio')
stopifnot(nrow(zfish_paths) >= 2)
query <- unique(as.character(getXrefList(zfish_paths$id[[1]], 'L')))
universe <- unique(unlist(lapply(zfish_paths$id[seq_len(min(20, nrow(zfish_paths)))], function(id) {
  as.character(getXrefList(id, 'L'))
})))
query <- query[query %in% universe]
stopifnot(length(query) >= 3, length(universe) > length(query))

ora <- enrichWP(query, organism = 'Danio rerio', universe = universe,
                minGSSize = 3, maxGSSize = 500, pvalueCutoff = 1.1, qvalueCutoff = 1.1)
stopifnot(nrow(as.data.frame(ora)) > 0)
cat('WikiPathways standalone API smoke passed:', length(organisms), 'organisms;',
    nrow(zfish_paths), 'Danio rerio pathways;', nrow(as.data.frame(ora)), 'ORA terms\n')
