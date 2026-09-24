# Purpose: MetaboAnalystR API-path ORA with a metabolome-filter background. SENDS the compound list to
#          https://www.xialab.ca/api/pathwayora (see "Undisclosed remote call" in SKILL.md); use only
#          if that is acceptable, otherwise use scripts/local_ora.R.
# Inputs:  compounds file (names/IDs, one per line); reference file (one KEGG ID per line); optional
#          organism (hsa), ID type (name), output CSV (ora_api.csv).
# Output:  the CSV (Raw p, FDR, Impact, Hits, Total per pathway), or a printed failure reason and
#          exit status 2: the server has been observed to reject the FILTERED request
#          (CalculateOraScore returns 0).
# Usage:   Rscript scripts/ora_api.R compounds.txt reference_metabolome.txt [hsa] [name] [ora_api.csv]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('usage: Rscript ora_api.R compounds.txt reference_metabolome.txt [organism] [idtype] [out]')
organism <- if (length(args) >= 3) args[3] else 'hsa'
id_type  <- if (length(args) >= 4) args[4] else 'name'
out_file <- if (length(args) >= 5) args[5] else 'ora_api.csv'

current.msg <- character(0); err.vec <- character(0)  # required -- see Version Compatibility
suppressMessages(library(MetaboAnalystR))

# Same mapping steps as scripts/map_compounds.R (the API call needs the mSet, not just the IDs)
mSet <- InitDataObjects('conc', 'pathora', FALSE)
mSet <- SetOrganism(mSet, organism)
compounds <- readLines(args[1]); compounds <- trimws(compounds[nzchar(trimws(compounds))])
mSet <- Setup.MapData(mSet, compounds)
mSet <- CrossReferencing(mSet, id_type)
mSet <- CreateMappingResultTable(mSet)

mSet <- SetKEGG.PathLib(mSet, organism, 'current')

# SetMetabolomeFilter(mSet, TRUE) alone does NOT restrict the background --
# Setup.KEGGReferenceMetabolome() must be called first to load the reference file into
# mSet$dataSet$metabo.filter.kegg, or the filter silently has no effect. FALSE uses the
# whole library (all of KEGG) -- the inflated default that manufactures false positives.
mSet <- SetMetabolomeFilter(mSet, TRUE)
mSet <- Setup.KEGGReferenceMetabolome(mSet, args[2])  # one KEGG ID per line

mSet <- CalculateOraScore(mSet, 'rbc', 'hyperg') # node-importance 'rbc'|'dgr'; test 'hyperg'|'fisher'
# Checked on MetaboAnalystR 4.3.0: the server has been observed to reject the FILTERED
# request outright (CalculateOraScore returns 0; current.msg == "Failed to connect to
# the API Server!"), even with a correctly-matched reference file, while the unfiltered
# (FALSE) call to the same endpoint succeeds. If this happens, use Local-Only ORA.
if (is.numeric(mSet)) {
  cat('ORA failed:', paste(current.msg, collapse = ' | '), '\n')
  quit(status = 2)
} else {
  ora <- as.data.frame(mSet$analSet$ora.mat)   # columns include Raw p, FDR, Impact, Hits, Total
  write.csv(ora, out_file)
  print(head(ora, 10))
}
