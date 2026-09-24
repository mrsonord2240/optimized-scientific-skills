# TMT / iTRAQ isobaric labeling

Read this when the data are TMT/iTRAQ: reporter extraction from mzML, CoA impurity correction with the orientation guard, and the multi-plex MSstatsTMT route with the reference-channel bridge.

### TMT/iTRAQ Isobaric Labeling
Reporter extraction is a spectra-level step, not a text-matrix read. Within a single plex the channels are co-isolated/co-fragmented in the same MS2 event, so relative ratios are stable; but MULTI-batch TMT CANNOT be compared across plexes without an IRS bridge (a pooled reference channel in every plex; Plubell 2017). Route to proteomics/quantification for the mechanics.
Run `scripts/tmt_impurity_correct.R tmt.mzML tmt10_coa.csv tmt_reporters_corrected.csv`. It builds the impurity matrix by channel name (never `makeImpuritiesMatrix(filename = )` for TMT10/TMTpro), stops on a transposed or fraction-scaled CoA, and refuses negative corrected values.

Multi-plex TMT (the common case: two or more plexes) -- do NOT concatenate plexes directly. MSstatsTMT applies the reference-channel (IRS-style) bridge during summarization. MaxQuant route, checked on MSstatsTMT 2.14.2 against its bundled 5-plex `evidence` / `proteinGroups` / `annotation.mq`.
Run `scripts/msstatstmt_multiplex.R evidence.txt proteinGroups.txt annotation.csv msstatstmt_results.csv`. The annotation needs a `Norm` reference channel in every plex; contrasts are built from the condition levels and re-adjusted across rows (`adj.pvalue.global`).
