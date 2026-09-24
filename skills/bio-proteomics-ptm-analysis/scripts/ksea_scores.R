#!/usr/bin/env Rscript
# Kinase activity from protein-adjusted site fold-changes with KSEAapp::KSEA.Scores.
#
# Builds the six-column PX table (Protein, Gene, Peptide, Residue.Both, p, FC) from the adjusted site table,
# guards the three known failure shapes (non-finite log2FC, empty PX, prior coverage < 2), then scores.
#
# Inputs (key=value):
#   adjusted=<csv>         adjusted_sites.csv written by msstatsptm_labelfree.R (columns Protein, log2FC, adj.pvalue);
#                          its contrast must be 'Treatment vs Control' (FC = 2^log2FC); if the Label reads
#                          'Control vs Treatment' negate log2FC first or every kinase's sign inverts
#   proteinGroups=<txt>    MaxQuant proteinGroups.txt of the global run (Protein IDs, Gene names)
#   prior=<csv>            FULL PhosphoSitePlus + NetworKIN table (KINASE, SUB_GENE, SUB_MOD_RSD, Source,
#                          networkin_score; github.com/casecpb/KSEA), not the abbreviated data(KSData)
#   out=<csv>              kinase scores (default kinase_scores.csv)
# Usage: Rscript scripts/ksea_scores.R adjusted=adjusted_sites.csv proteinGroups=proteinGroups_global.txt \
#          prior='PSP&NetworKIN_Kinase_Substrate_Dataset.csv' out=kinase_scores.csv
opt <- list(adjusted = 'adjusted_sites.csv', proteinGroups = 'proteinGroups_global.txt',
            prior = 'PSP&NetworKIN_Kinase_Substrate_Dataset.csv', out = 'kinase_scores.csv')
for (a in strsplit(commandArgs(trailingOnly = TRUE), '=', fixed = TRUE)) opt[[a[1]]] <- paste(a[-1], collapse = '=')
rd <- function(f) read.table(f, sep = '\t', header = TRUE, quote = '')
adjusted <- read.csv(opt$adjusted)

library(KSEAapp)
KSData <- read.csv(opt$prior)   # user-supplied prior
pg <- rd(opt$proteinGroups)
gene_symbol <- setNames(sub(';.*', '', pg$Gene.names), sub(';.*', '', pg$Protein.IDs))

# log2FC = +/-Inf (site measured in one condition only) gives FC = 0 or Inf, and one such row
# turns every kinase z-score NaN; filter on log2FC, not on the linear FC (0 is finite)
ks <- adjusted[is.finite(adjusted$log2FC), ]
if (nrow(ks) == 0)
  stop('No site has a finite log2FC -- every row of the adjusted table is +/-Inf or NA, i.e. every ',
       'site was quantified in one condition only. KSEA scores fold changes, not presence/absence; ',
       'report those sites as detected-in-one-condition instead.')
PX <- data.frame(
  Protein = sub('_[STY][0-9]+$', '', ks$Protein),
  Gene = gene_symbol[sub('_[STY][0-9]+$', '', ks$Protein)],   # HUGO symbol; merge key with SUB_GENE
  Peptide = rep('NULL', nrow(ks)),   # rep(), not the bare literal: a length-1 value against
                                     # zero-length columns is 'differing number of rows: 0, 1'
  Residue.Both = sub('^.*_', '', ks$Protein),                  # e.g. S473; merge key with SUB_MOD_RSD
  p = ks$adj.pvalue,
  FC = 2^ks$log2FC)   # Treatment/Control because the contrast is 'Treatment vs Control'
PX <- PX[!is.na(PX$Gene), ]
if (nrow(PX) == 0)
  stop('PX is empty after dropping sites with no gene symbol: proteinGroups_global.txt carried no ',
       '`Gene names` for any tested protein. Check the search FASTA had gene annotation, or map the ',
       'accessions to HUGO symbols yourself before building PX.')

# Check PRIOR COVERAGE before calling. KSEA.Scores merges on SUB_GENE + SUB_MOD_RSD and then
# aggregates, so a prior overlapping this site list in 0 or 1 place dies inside the package with
# `no rows to aggregate` -- a message that says nothing about coverage. Curated priors cover only a
# small fraction of any real site list, so this is the common failure, not an exotic one.
prior <- KSData[grep('PhosphoSitePlus', KSData$Source), ]   # the subset NetworKIN = FALSE will use
covered <- sum(paste(PX$Gene, PX$Residue.Both) %in% paste(prior$SUB_GENE, prior$SUB_MOD_RSD))
cat('sites in PX:', nrow(PX), '| covered by the prior:', covered, '\n')
if (covered < 2)
  stop('The kinase-substrate prior covers ', covered, ' of ', nrow(PX), ' sites, so KSEA has ',
       'nothing to score. Check that KSData is the FULL PhosphoSitePlus+NetworKIN table (not the ',
       'abbreviated data(KSData)), that Gene holds HUGO symbols matching SUB_GENE, and that ',
       'Residue.Both is formatted like SUB_MOD_RSD (S473, not pS473 or Ser473).')

# NetworKIN = FALSE: PhosphoSitePlus-curated pairs only; TRUE adds predictions above NetworKIN.cutoff
# (and then coverage should be counted against that subset instead).
kinase_scores <- KSEA.Scores(KSData, PX, NetworKIN = FALSE, NetworKIN.cutoff = 3)
kinase_scores <- kinase_scores[order(kinase_scores$z.score), c('Kinase.Gene', 'm', 'z.score', 'FDR')]
print(kinase_scores)
write.csv(kinase_scores, opt$out, row.names = FALSE)
