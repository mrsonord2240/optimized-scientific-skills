# Purpose: map identified compound names/IDs to KEGG compounds with MetaboAnalystR and print the
#          mapping table, so coverage is inspected before any enrichment p-value is trusted.
# Inputs:  compounds file (one name/ID per line); optional organism (default hsa), ID type
#          (name|hmdb|kegg|pubchem, default name), output file (default kegg_ids.txt).
# Output:  the output file, one mapped KEGG compound ID per line.
# Usage:   Rscript scripts/map_compounds.R compounds.txt [hsa] [name] [kegg_ids.txt]
# Note:    the first run downloads MetaboAnalyst's generic reference libraries into the working
#          directory (see Version Compatibility in SKILL.md); no compound list is sent.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop('usage: Rscript map_compounds.R compounds.txt [organism] [idtype] [out]')
compounds_file <- args[1]
organism <- if (length(args) >= 2) args[2] else 'hsa'
id_type  <- if (length(args) >= 3) args[3] else 'name'
out_file <- if (length(args) >= 4) args[4] else 'kegg_ids.txt'

current.msg <- character(0); err.vec <- character(0)  # required -- see Version Compatibility
library(MetaboAnalystR)

# 'pathora' = pathway ORA; 'conc' = concentration-style input
mSet <- InitDataObjects('conc', 'pathora', FALSE)
mSet <- SetOrganism(mSet, organism)

# Confidently identified compounds (MSI level 1-2); names, HMDB, or KEGG IDs
compounds <- readLines(compounds_file)
compounds <- trimws(compounds[nzchar(trimws(compounds))])
mSet <- Setup.MapData(mSet, compounds)
mSet <- CrossReferencing(mSet, id_type)         # 'name' | 'hmdb' | 'kegg' | 'pubchem'
mSet <- CreateMappingResultTable(mSet)          # inspect mapping coverage before trusting any p-value
kegg_ids <- mSet$dataSet$map.table[, 'KEGG']
kegg_ids <- kegg_ids[!is.na(kegg_ids) & nzchar(kegg_ids)]

cat(sprintf('Mapped %d of %d compounds to KEGG IDs\n', length(kegg_ids), length(compounds)))
writeLines(kegg_ids, out_file)
