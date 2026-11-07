# Reference Mapping Failure Modes

Reference-mapping failure modes: scANVI/kNN forcing labels, softmax versus transfer uncertainty, label leakage, reference composition bias, query QC artifacts, gene-set mismatch, integration scores versus labels, zero-shot foundation models; reconciliation when methods disagree.

## Per-Method Failure Modes

### scANVI / kNN -- forcing the query onto reference labels
- **Trigger:** Query contains a population absent from the reference (novel type, disease state, perturbed program).
- **Mechanism:** The classifier/kNN assigns every query cell to its nearest reference label; there is no "none of the above" unless added.
- **Symptom:** A coherent novel cluster split across 2-3 reference labels, each with *high* probability; the UMAP looks "well integrated."
- **Fix:** Gate on transfer uncertainty / OOD distance (above); this catches types far from the reference but not a missing type adjacent to a reference type (measured: monocytes removed -> 95.9% called DC, 1.6% gated). Run the per-label marker check. Inspect query-only clusters for marker genes independent of transferred labels.

### Softmax overconfidence vs label-transfer uncertainty conflated
- **Trigger:** Reporting "confidence" as the `predict(soft=True)` max.
- **Mechanism:** The softmax measures *which* reference label conditional on belonging; it is normalized away from distance and cannot say "far from everything." A cell can be 0.99 "T cell" and be a hepatocyte.
- **Symptom:** Pipeline filters on softmax >= 0.5 and still passes OOD cells.
- **Fix:** Threshold the weighted-kNN uncertainty (HLCA 0.2) or a Mahalanobis/ensemble OOD signal for the "Unknown" decision; use the softmax only to disambiguate among in-distribution labels.

### scANVI -- semi-supervised label leakage / latent carving
- **Trigger:** scANVI reference where labels strongly drive latent geometry; trajectory or novel-state query.
- **Mechanism:** The classifier head back-propagates label structure into the latent, carving it to separate reference types; query cells are pulled toward that structure even when their biology lies between/outside it.
- **Symptom:** A query continuum (differentiation trajectory) collapses onto discrete reference clusters; intermediate states vanish.
- **Fix:** For trajectory/novel-state queries prefer *unsupervised* scVI surgery, annotate the query independently, and cross-check against the scVI latent.

### Reference composition / missing-biology bias
- **Trigger:** Reference from healthy/limited donors; query from disease, different age, ancestry, or tissue region.
- **Mechanism:** The reference manifold is the entire hypothesis space; off-manifold cells are projected onto the nearest on-manifold point and rare reference populations act as attractors.
- **Symptom:** Disease-specific states labeled as the closest healthy type; ancestry/age effects read as "batch."
- **Fix:** Audit reference composition before mapping; prefer references covering the query's expected biology; treat mapping as hypothesis generation and validate query findings de novo; consider extending the reference (treeArches).

### Query QC artifacts laundered into confident labels
- **Trigger:** Query not QC'd to the reference's standard (empty droplets, ambient RNA, doublets, high-MT).
- **Mechanism:** A doublet sits between two reference types and maps to a spurious "intermediate"; ambient RNA shifts profiles toward the dominant type.
- **Symptom:** Artifactual "transitional" populations; doublet clusters labeled as rare real types.
- **Fix:** Run full query QC *before* mapping (single-cell/doublet-detection, ambient correction, MT/count filters matched to the reference). Mapping does not clean data.

### Feature-space / gene-set mismatch
- **Trigger:** Query missing reference HVGs; different gene annotation/version; `prepare_query_anndata` skipped.
- **Mechanism:** The encoder expects the exact reference gene order; missing genes are zero-padded and reordering silently corrupts the input.
- **Symptom:** Garbage latent, everything maps to one blob, or a silent accuracy cliff (no error raised).
- **Fix:** Always `prepare_query_anndata(query, reference_model)`; verify the shared-gene fraction; too few shared HVGs is a hard stop.

### Good integration metrics, wrong labels
- **Trigger:** Judging mapping by scIB integration scores alone.
- **Mechanism:** Integration metrics reward query/reference mixing; mixing OOD cells into the wrong neighborhood *raises* the batch-removal score, and bio-conservation uses reference labels (circular for the query).
- **Symptom:** Top scIB total score with biologically wrong annotation.
- **Fix:** Integration metrics validate the *embedding*, not labels. Evaluate transfer on held-out labeled query cells (per-type F1, especially rare types), OOD detection on spiked-in unseen types, and marker-gene sanity checks.

### Zero-shot foundation-model embedding as a mapper
- **Trigger:** Using scGPT/Geneformer zero-shot embeddings for clustering/transfer expecting SOTA.
- **Mechanism:** The masked-gene pretraining objective does not guarantee a label- or batch-aware embedding; zero-shot embeddings are not batch-corrected.
- **Symptom:** Worse AvgBio/integration than scVI or even HVG-PCA + Harmony.
- **Fix:** Fine-tune with task labels, or use an established mapper; reserve foundation models for cross-modality/species/data-scarce cases (Kedzierska 2025).

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| scANVI label confident but kNN uncertainty high | OOD cell forced onto nearest label | Trust the uncertainty; set Unknown and inspect markers |
| Symphony Mahalanobis flags OOD but scANVI does not | scANVI latent carved to absorb the cell | Prefer the distance-based flag; novel biology likely |
| popV members disagree | Genuine ambiguity or granularity mismatch | Route to manual review; report the disagreement, do not force a leaf |
| High scIB score, poor per-type F1 on held-out labels | Embedding mixes well but labels wrong | Believe the F1; integration score is not a label metric |
