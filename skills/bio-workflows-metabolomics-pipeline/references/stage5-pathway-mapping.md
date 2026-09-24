# Stage 5 -- Pathway Mapping (the background is the null)

**Goal:** Interpret the differential result in pathway context without laundering annotation uncertainty into confident biology.

**Approach:** Two disjoint entry points. Confidently identified compounds -> ORA/MSEA with an assay-coverage background (NOT all of KEGG). Raw m/z with no IDs -> mummichog/PSEA whose permutation null is sampled from the FULL feature table. Either way, report mapping coverage and the MSI levels of the driving compounds; downgrade claims to "consistent with perturbation." Full method choice and background construction in metabolomics/pathway-mapping.

```r
# MSI/Schymanski gate (commitment #2): only Level 1-2 (1, 2a, 2b) may enter identified-ORA as an
# "identification" -- Level 3-5 is a database-name hypothesis, not a name (see
# metabolite-annotation). This filter is the one place the principle is enforced as code, not
# just stated in a QC-checkpoint row.
# annotated: Stage 3's per-feature output (feature_id, name, msi_level -- msi_level values like
# 1, '2a', '2b', 3, 4, 5 per metabolite-annotation's assign_level()), joined to feature_id.
identified_compounds <- unique(annotated$name[grepl('^[12]', as.character(annotated$msi_level))])
stopifnot(length(identified_compounds) > 0)   # all Level 3-5? use Path B (mummichog) instead, not a forced ORA

current.msg <- character(0); err.vec <- character(0)  # required outside the Shiny app; see pathway-mapping
library(MetaboAnalystR)
# Path A: identified compounds (MSI level 1-2) -> ORA
mSet <- InitDataObjects('conc', 'pathora', FALSE)
mSet <- SetOrganism(mSet, 'hsa')
mSet <- Setup.MapData(mSet, identified_compounds)
mSet <- CrossReferencing(mSet, 'name')
mSet <- CreateMappingResultTable(mSet)              # inspect coverage before trusting any p-value
mSet <- SetKEGG.PathLib(mSet, 'hsa', 'current')
mSet <- SetMetabolomeFilter(mSet, TRUE)             # TRUE alone does not restrict the background --
mSet <- Setup.KEGGReferenceMetabolome(mSet, 'reference_metabolome.txt')  # this call does; file = one KEGG compound ID per line, names fail (see pathway-mapping)
mSet <- CalculateOraScore(mSet, 'rbc', 'hyperg')
# Sends the mapped compound list to xialab.ca and can reject a filtered request outright on this
# version; if so, use pathway-mapping's Local-Only ORA (KEGGREST + local phyper, no remote call).

# Path B: no IDs -> mummichog on the FULL peak table (m/z + p-value + t-score)
# mSet <- InitDataObjects('mass_all', 'mummichog', FALSE)
# mSet <- UpdateInstrumentParameters(mSet, 5.0, 'negative')   # ppm + ionization mode are mandatory
# mSet <- Read.PeakListData(mSet, 'peaks_full.txt')           # ENTIRE table, not significant-only
# mSet <- PerformPSEA(mSet, 'hsa_mfn', 'current', permNum = 1000)
```
