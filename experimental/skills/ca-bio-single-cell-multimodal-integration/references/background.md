# Background (not needed to run an analysis)

Anchor structure: horizontal (same modality, different cells; anchor shared features), vertical (paired, same cell; anchor shared cells), diagonal (different modalities in different cells; correspondence inferred from prior knowledge), mosaic (partially observed grid of modalities by batches).

CITE-seq ADT background is a three-part mixture: ambient antibody captured in every droplet including empties; cell-intrinsic non-specific binding (Fc receptors, sticky dying cells) that does not appear in empties; spillover and index hopping between barcodes.

WNN weights reward local neighbour predictability, so a handful of high-variance or saturating ADT features can manufacture self-consistent neighbourhoods and be up-weighted despite carrying less biology.

Imputed modalities: MultiVI, StabMap and Cobolt impute the missing modality for unpaired cells; gene-activity scores from ATAC approximate RNA.

Model notes: WNN learns per-cell, per-modality weights from cross-modality neighbour prediction. totalVI is a conditional VAE with RNA NB/ZINB and each protein a two-component NB mixture. MultiVI is a joint VAE over RNA, ATAC and protein, mosaic-capable. MOFA+ is a linear Bayesian group factor analysis. mojitoo is CCA across precomputed reductions. GLUE is per-modality VAEs with a prior feature graph and adversarial alignment. The bridge is a multiome dictionary linking the query modality to the reference. StabMap projects through shortest paths over shared features.

Source test results: totalVI/MultiVI synthetic checks (RNA-only cells land nearer same-type paired cells, 0.27 against 0.56); bridge synthetic recovery 81% and 75%.

The source `examples/cite_seq_analysis.R` and `.py` became `scripts/cite_seq_wnn.R` and `scripts/cite_seq_muon.py`; the Multiome WNN, MOFA+ and GLUE snippets became `scripts/multiome_wnn.R`, `scripts/mofa_muon.py` and `scripts/glue_unpaired.py`. The source's prerequisites were `install.packages(c('Seurat','dsb'))`, `BiocManager::install('Signac')` and `pip install muon mudata scanpy anndata scvi-tools scglue`.
