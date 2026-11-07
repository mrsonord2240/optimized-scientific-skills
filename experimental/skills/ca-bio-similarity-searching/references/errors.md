## Per-Tool Failure Modes

### ECFP4 Tanimoto -- saturation in diverse libraries

**Trigger:** Library spans drug-like + natural products + peptides + metabolites.

**Mechanism:** A local circular fingerprint may not preserve the distinctions needed for a particular mixed-modality retrieval task.

**Symptom:** Known relevant neighbors are not enriched above background, or retrieval metrics and neighborhood stability are poor on target-relevant controls. A low mean pairwise similarity alone is not a universal saturation test.

**Fix:** Benchmark ECFP4 against alternatives such as MHFP6 or MAP4 using held-out analog recovery, scaffold-aware retrieval, or another task-aligned metric. Do not transfer raw-score thresholds between fingerprint families.

### Butina clustering -- O(N^2) memory blowup

**Trigger:** Library >100k molecules.

**Mechanism:** Butina requires upper-triangle distance matrix, ~5e9 floats for 100k compounds.

**Symptom:** OOM error or hours of CPU.

**Fix:** Use approximate clustering (HDBSCAN on UMAP-reduced fingerprints) or LSH-based clustering on MHFP6.

### MCS -- timeout versus a genuinely trivial common fragment

**Trigger:** Mol set with low overlap, large molecules, or many input mols.

**Mechanism:** MCS search can be combinatorial for difficult mappings, but divergent molecules can also converge quickly on a genuinely small common fragment.

**Symptom:** Returns a small partial MCS or empty result. Distinguish the two cases with `result.canceled`.

**Fix:** If `result.canceled=True`, reduce input count or raise `timeout` only if the cost is acceptable. If `result.canceled=False` and the MCS is small, do not expect a longer timeout to help: pre-cluster by Tanimoto and run MCS within coherent clusters.

### Tanimoto = 1.0 != same molecule

**Trigger:** Comparing stereoisomers or distinct molecules with identical fingerprint bits.

**Mechanism:** RDKit Morgan fingerprints use `includeChirality=False` unless requested, so enantiomers and some diastereomers can be indistinguishable. Folded hashed fingerprints can also collide, especially as molecules grow or bit length shrinks.

**Symptom:** Two structurally different molecules report Tanimoto 1.0.

**Fix:** For stereochemistry-sensitive similarity, build both query and library fingerprints with `GetMorganGenerator(..., includeChirality=True)`. For exact identity, compare standardized canonical SMILES or InChIKey, not any fingerprint. An unhashed sparse fingerprint can help diagnose hash collisions but does not rescue a chirality-blind representation.

### Similarity threshold transfer fails

**Trigger:** Threshold tuned on ECFP4 applied to RDKit FP, AtomPair, or MACCS.

**Mechanism:** Bit-density and fragment-resolution differ; Tanimoto distributions shift.

**Symptom:** "Similar" set is much larger or smaller than expected.

**Fix:** Re-tune the threshold per fingerprint and dataset. Compute its pairwise-similarity distribution and choose a retained-pair percentile that matches the review capacity (for example p99.5); sample pairs for large libraries. AtomPair ~0.55, MACCS ~0.85, ECFP4 ~0.7, and FCFP4 ~0.6 are repository starting heuristics, not universal equivalents.

## Reconciliation: Cliffs Across Methods

If a pair flags as an activity cliff under one representation but not another, treat that as representation sensitivity. Inspect atom mappings, fingerprint environments, assay uncertainty, and the exact structural change; disagreement alone does not establish which substituent caused the activity difference.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `BulkTanimotoSimilarity` output is treated as bit counts | The API returns similarity values, normally floats | Keep the returned values as similarities; inspect input vector types if the output is unexpected |
| Reported similarity > 1 | Custom formula, malformed data, negative features, or an incorrectly normalized external implementation | Verify the coefficient definition and inputs; standard nonnegative RDKit Tanimoto and Tversky similarities are bounded by 1 |
| Cluster centroids change when input order changes | Taylor-Butina tie handling and assignment are order-sensitive | Standardize and sort inputs by a stable identifier before clustering; record `reordering` and the input order |
| MaxMinPicker returns first N inputs | All-zero initial similarity matrix | Seed picker explicitly: `picker.LazyBitVectorPick(fps, n_lib, n_pick, seed=42)` |
| Activity cliff "false positives" | Chirality-blind Morgan fingerprints merge stereoisomers; unstripped salts can create duplicate-like pairs; hash collisions are secondary | Standardize and strip salts first; use `includeChirality=True` when stereochemistry matters; compare standardized canonical SMILES or InChIKey for identity |
| Diverse subset has duplicates | Standardization not applied | Canonicalize via `chemoinformatics/molecular-standardization` first |
| Tanimoto incompatible with neural fingerprint | Continuous-valued fingerprint | Use cosine or sklearn `cdist` with `'cosine'` metric |
