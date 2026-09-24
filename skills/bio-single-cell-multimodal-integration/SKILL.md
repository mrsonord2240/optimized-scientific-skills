---
name: bio-single-cell-multimodal-integration

description: Integrate multimodal single-cell data (CITE-seq RNA+protein, 10x Multiome RNA+ATAC, unpaired/diagonal RNA+ATAC) and choose the right joint method. Use when classifying an integration task by anchor structure (paired vs unpaired), denoising CITE-seq ADT background before joint embedding, picking between WNN, totalVI, MultiVI, MOFA+, GLUE, or Seurat v5 bridge integration, or diagnosing why a modality dominates a joint clustering.

tool_type: mixed

primary_tool: Seurat

license: MIT

author: GPTomics

---



## Version Compatibility

Reference examples tested with: scanpy 1.10+, Seurat 5.0+, anndata 0.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Muon WNN is stored on a `MuData` object, not an `AnnData` object: after
`mu.pp.neighbors(..., key_added='wnn')`, run
`mu.tl.umap(mdata, neighbors_key='wnn')`. Do not pass Muon's multimodal WNN
metadata to `scanpy.tl.umap`, whose AnnData neighbor contract is different.

## Prerequisites

```r
install.packages(c('Seurat', 'dsb'))
BiocManager::install('Signac')          # Multiome ATAC
```

```bash
pip install muon mudata scanpy anndata scvi-tools
pip install scglue                       # unpaired/diagonal integration; no Windows install (pybedtools -> pysam has no wheel, pip fails building it): use WSL/Linux/macOS, or skip if GLUE is not needed
```

# Multimodal Integration

**"Jointly analyze my CITE-seq / Multiome / unpaired multi-omic data"** -> Classify the task by anchor structure, denoise each modality in its native pipeline, then build one joint representation.
- R: `Seurat::FindMultiModalNeighbors()` (WNN), `Signac` (ATAC LSI), `dsb::DSBNormalizeProtein()` (ADT denoising), `PrepareBridgeReference()` (v5 bridge)
- Python: `muon`/`mudata` (MuData container), `scvi.model.TOTALVI` / `MULTIVI`, `MOFA2`/`muon.tl.mofa`, `scglue` (diagonal)

## Governing Principle

Classify the integration task by its anchor structure FIRST, because the anchor decides which algorithm class is even applicable (Argelaguet 2021).
- Horizontal: same modality, different cells; anchor = shared features (batch correction, not this skill).
- Vertical (paired, same cell): multiple modalities measured in the same cells; anchor = shared cells. CITE-seq, 10x Multiome.
- Diagonal (unpaired): different modalities in different cells, no shared cells and no shared features; correspondence is inferred from prior knowledge. Independent scRNA + scATAC.
- Mosaic: partially observed grid of (modalities x batches); some blocks present, some missing.

Paired vs unpaired is the master fork: paired correspondence is known a priori (WNN, totalVI, MultiVI, MOFA+), unpaired/diagonal correspondence must be inferred (GLUE, Seurat v5 bridge), mosaic mixes both (MultiVI, StabMap). Two separately-paired datasets that share only one modality (for example a 10x Multiome and a CITE-seq experiment sharing only RNA) are a mosaic problem: anchor on the shared RNA and impute or bridge the modality-specific blocks with StabMap, MultiVI, or Seurat v5 bridge integration rather than forcing a single WNN.

CITE-seq ADT background is a three-part mixture, not one "ambient" term: (1) ambient antibody captured in every droplet including empties, (2) cell-intrinsic non-specific binding (Fc receptors, sticky dying cells) that does NOT appear in empties, (3) spillover/index hopping between barcodes. Denoise ADT (DSB or totalVI's built-in background mixture) BEFORE any joint embedding; raw or CLR-only ADT carries this background into the joint graph.

WNN can be dominated by the noisier modality: weights reward local neighbor predictability, and a handful of high-variance or saturating ADT features can manufacture self-consistent neighborhoods and get up-weighted despite carrying less biology. Report the per-cell weight distribution and check whether clustering survives down-weighting the suspect modality.

Imputed modalities are inferences, not measurements: MultiVI/StabMap/Cobolt impute the missing modality for unpaired cells, and gene-activity scores from ATAC approximate RNA. Differential expression or marker calls on imputed values are model-dependent and must be flagged as such.

## Classify the Task: Anchor Structure -> Method Class

| Anchor structure | What is shared | Example assay | Method class |
|---|---|---|---|
| Vertical / paired | Same cells | CITE-seq, 10x Multiome | WNN, totalVI, MultiVI(paired), MOFA+, mojitoo |
| Diagonal / unpaired | Nothing (prior graph) | Independent scRNA + scATAC | GLUE, Seurat v5 bridge, LIGER |
| Mosaic | Some modalities only | Batch A RNA+ATAC, batch B RNA | MultiVI, StabMap, Cobolt, totalVI(partial) |

When methods compete, verify the current best-practice default against the installed tool docs before committing; the field moves and defaults drift across minor versions.

## Method Decision Table (Paired CITE-seq / Multiome)

| Method | Model / assumption | Use when | Fails when |
|---|---|---|---|
| WNN (Seurat; `references/cite-seq-dsb-wnn.md`, `references/multiome-mofa.md`) | Per-cell, per-modality weights from cross-modality neighbor prediction; one weighted graph | Fast joint embedding/clustering of one well-normalized paired dataset | Protein background not removed upstream; noisy/saturating modality dominates; not for unpaired/mosaic |
| totalVI (scvi-tools; `references/scvi-totalvi-multivi.md`) | Conditional VAE; RNA NB/ZINB, each protein a 2-component NB mixture (background+foreground) | Need denoised protein, principled DE, batch integration, merging different antibody panels | Tiny datasets (VAE overfits); no GPU and very large data; protein-specific background structure not captured by one per-cell factor |
| MultiVI (scvi-tools; `references/scvi-totalvi-multivi.md`) | Single joint VAE over RNA+ATAC(+protein); mosaic-capable, imputes missing modality | Paired+unpaired RNA/ATAC mixed (mosaic); want generative DE/DA | "batch" key is the modality indicator, not sequencing batch; imputed modalities treated as measured |
| MOFA+ (MOFA2; `references/multiome-mofa.md`) | Linear Bayesian group factor analysis; sparse factors, per-modality variance explained | Interpreting shared vs modality-specific axes of variation (exploratory/explanatory) | Used for clustering/denoising; likelihood mismatched to data; expecting batch correction within a view |
| mojitoo | CCA across precomputed per-modality reductions; fast, parameter-free | Quick paired joint reduction from existing PCA/LSI slots | No knob to down-weight a noisy modality; bounded by input reductions; paired only |

## Method Decision Table (Unpaired / Diagonal / Mosaic)

| Method | Model / assumption | Use when | Fails when |
|---|---|---|---|
| GLUE (scglue; `references/unpaired-glue-bridge.md`) | Per-modality VAEs + prior feature graph (peak-near-gene); adversarial cell alignment | Unpaired diagonal scRNA + scATAC; want regulatory inference as a byproduct | Genome-build/coordinate mismatch yields an empty guidance graph and garbage alignment; adversarial over-mixing of distinct states |
| Seurat v5 bridge (`references/unpaired-glue-bridge.md`) | Multiome bridge dataset = dictionary linking query modality to reference modality | Mapping a query (scATAC) onto a reference built in another modality (scRNA) | Poor/batch-mismatched bridge propagates error; rare query-only populations mislabeled |
| StabMap | Mosaic topology from shared features; project all cells via shortest paths | Mosaic with informative unshared features that cannot be dropped | Unshared-feature chaining compounds error per hop |
| Cobolt / scMoMaT | Generative shared latent over joint + single-modality datasets | Mosaic where a generative latent is preferred over feature chaining | DE/marker calls made on imputed values |

## ADT Normalization: CLR vs DSB

| Method | What it does | Use when | Fails when |
|---|---|---|---|
| CLR (centered log-ratio) | Rescales compositionally; Seurat `NormalizeData(method="CLR", margin=2)` | Quick, no empty droplets available; small panels | Does NOT remove background; geometric-mean denominator distorted by saturating high-abundance ADTs |
| DSB (`references/cite-seq-dsb-wnn.md`) | Ambient correction from empty droplets + per-cell technical denoising via 2-component mixture + isotype controls | Raw/unfiltered matrix available (needs empty droplets); want background removed before embedding | No empty droplets retained; protein-specific non-specific binding (one per-cell factor under/over-corrects); no clearly bimodal proteins |

Seurat's CLR margin is genuinely ambiguous across versions (margin=2 = per-feature is the WNN-tutorial recommendation for large panels); verify with `?NormalizeData` on the installed version.

## Reference Files

Read the file for the method you are running; `SKILL.md` above decides which method applies.

| File | Read when |
|---|---|
| `references/cite-seq-dsb-wnn.md` | CITE-seq: denoising ADT with DSB from the raw matrix, then WNN joint clustering (Seurat). Runnable end to end: `examples/cite_seq_analysis.R` |
| `examples/cite_seq_analysis.py` | CLR-only Python/Muon CITE-seq WNN fallback. It creates modality-local Scanpy graphs, then uses Muon's WNN and MuData-aware UMAP path. |
| `references/scvi-totalvi-multivi.md` | Training totalVI (CITE-seq denoising + DE) or MultiVI (mosaic, RNA+ATAC partially observed) with scvi-tools. Runs `scripts/totalvi_cite_seq.py`, `scripts/multivi_mosaic.py` |
| `references/multiome-mofa.md` | 10x Multiome RNA + ATAC WNN (Signac LSI), or MOFA+ shared/specific factors |
| `references/unpaired-glue-bridge.md` | Unpaired scRNA + scATAC: GLUE, or Seurat v5 bridge integration through a multiome bridge. Runs `scripts/seurat_bridge_integration.R` |

## MuData Housekeeping

After per-modality QC, modalities hold different cell sets; `muon.pp.intersect_obs(mdata)` before any paired analysis. Editing a modality-local `mdata.mod['rna'].obs` needs `mdata.update()` to propagate to the global `mdata.obs`. R round-trips (MuDataSeurat, zellkonverter) are lossy; plan to stay in one ecosystem.

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| WNN clustering driven entirely by ADT | A few saturating high-variance proteins dominate the neighbor graph | Report per-cell weight distribution; down-weight or denoise ADT (DSB); re-check clustering stability |
| "Background" smear in every ADT cluster | Ran WNN/CLR without empty-droplet denoising | Run DSB (needs raw/unfiltered matrix) or totalVI before joint embedding |
| DSB errors / nonsense output | Passed a filtered cell matrix only (no empty droplets) | Supply `empty_drop_matrix` from the raw/unfiltered matrix |
| Spurious batch structure after merging multiome | Per-dataset peak sets, not a unified set | Re-quantify all cells against one common peak set |
| GLUE produces a blob / no alignment | Guidance graph near-empty from genome-build/coordinate mismatch | Align RNA gene coords and ATAC peaks to the same build before building the graph |
| RNA and protein disagree for a marker | Often real post-transcriptional biology (stability, trafficking, lag), not an artifact | Do not "correct away"; treat single-gene discordance as informative |
| MultiVI batch effects persist | The `batch_key` was set to the modality indicator, not sequencing batch | Add a separate covariate for the real batch |
| DE on a modality looks too clean | Computed on imputed/gene-activity values, not measurements | Flag imputed-modality DE as model-dependent; validate against a measured modality |

## Related Skills

- single-cell/scatac-analysis - ATAC QC, TF-IDF/LSI, gene-activity caveats for the Multiome ATAC half
- single-cell/preprocessing - per-modality RNA QC and normalization before integration
- single-cell/clustering - clustering and UMAP on the joint graph
- single-cell/batch-integration - horizontal (same-modality, cross-sample) correction
- single-cell/markers-annotation - marker-based interpretation of joint clusters
- atac-seq/motif-deviation - chromVAR TF activity on the Multiome ATAC modality
- pathway-analysis/go-enrichment - functional interpretation of modality-specific factors

## References

Argelaguet R, Cuomo ASE, Stegle O, Marioni JC. Computational principles and challenges in single-cell data integration. Nat Biotechnol 39(10):1202-1215 (2021).
Stoeckius M, Hafemeister C, Stephenson W, et al. Simultaneous epitope and transcriptome measurement in single cells (CITE-seq). Nat Methods 14:865-868 (2017).
Mulè MP, Martins AJ, Tsang JS. Normalizing and denoising protein expression data from droplet-based single-cell profiling (DSB). Nat Commun 13:2099 (2022).
Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data (WNN). Cell 184(13):3573-3587 (2021).
Gayoso A, Steier Z, Lopez R, et al. Joint probabilistic modeling of single-cell multi-omic data with totalVI. Nat Methods 18:272-282 (2021).
Ashuach T, Gabitto MI, Koodli RV, et al. MultiVI: deep generative model for the integration of multimodal data. Nat Methods 20(8):1222-1231 (2023).
Argelaguet R, Arnol D, Bredikhin D, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome Biol 21:111 (2020).
Cao Z-J, Gao G. Multi-omics single-cell data integration and regulatory inference with graph-linked unified embedding (GLUE). Nat Biotechnol 40(10):1458-1466 (2022).
Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis (Seurat v5 bridge). Nat Biotechnol 42:293-304 (2024).
Bredikhin D, Kats I, Stegle O. MUON: multimodal omics analysis framework. Genome Biol 23:42 (2022).
Ghazanfar S, Guibentif C, Marioni JC. Stabilized mosaic single-cell data integration using unshared features (StabMap). Nat Biotechnol 42(2):284-292 (2024).
Yin Y, et al. Characterization and decontamination of background noise in droplet-based single-cell protein expression data with DecontPro. Nucleic Acids Res 52(1):e4 (2024).
