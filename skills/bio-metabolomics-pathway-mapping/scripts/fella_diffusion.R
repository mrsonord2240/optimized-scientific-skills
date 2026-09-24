# Purpose: FELLA network-diffusion enrichment: the enzymes, reactions and modules that link the
#          affected metabolites (mechanism), not just a ranked pathway list.
# Inputs:  KEGG compound IDs file (one per line; KEGG IDs only); optional database dir (fella_hsa,
#          built once and reused if it exists); optional output CSV (fella_results.csv).
# Output:  the CSV from generateResultsTable(); unmapped compounds are printed (report them).
# Usage:   Rscript scripts/fella_diffusion.R cpd_ids.txt [fella_hsa] [fella_results.csv]
# Note:    first run builds the KEGG graph from the live KEGG API (about a minute) and the diffusion
#          matrices (about 10 minutes, ~630 MB on disk in the database dir).
suppressMessages(library(FELLA))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop('usage: Rscript fella_diffusion.R cpd_ids.txt [db_dir] [out.csv]')
db_dir   <- if (length(args) >= 2) args[2] else 'fella_hsa'
out_file <- if (length(args) >= 3) args[3] else 'fella_results.csv'

# Build once, reuse. buildGraphFromKEGGREST hits the live KEGG API (slow); cache the DB.
if (!dir.exists(db_dir)) {
  graph <- buildGraphFromKEGGREST(organism = 'hsa')
  buildDataFromGraph(keggdata.graph = graph, databaseDir = db_dir, internalDir = FALSE)
}
fella.data <- loadKEGGdata(databaseDir = db_dir, internalDir = FALSE)

cpd_ids <- readLines(args[1]); cpd_ids <- trimws(cpd_ids[nzchar(trimws(cpd_ids))]) # KEGG compound IDs only
analysis <- defineCompounds(compounds = cpd_ids, data = fella.data)
print(getExcluded(analysis))                       # compounds that did not map -- report this

# 'diffusion' is the recommended default; runHypergeom = plain ORA over the graph,
# runPagerank (lowercase r) = directed random walks. The method string is lowercase.
# approx = 'normality' (shown here) is analytic/deterministic, no seed needed. If using
# approx = 'simulation' instead, call set.seed() first -- it resamples niter times and is
# not reproducible run-to-run otherwise.
analysis <- runDiffusion(object = analysis, data = fella.data, approx = 'normality')
results <- generateResultsTable(object = analysis, data = fella.data, method = 'diffusion', threshold = 0.05)
write.csv(results, out_file)
print(head(results, 15))
