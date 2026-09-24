---
name: bio-crispr-screens-perturb-seq-analysis
category: Data Analysis
description: Analyzes single-cell pooled CRISPR screens (Perturb-seq, CROP-seq, Perturb-CITE-seq, ECCITE-seq, multiome) where each cell carries an sgRNA and a scRNA-seq / surface-protein / chromatin readout. Covers experimental design (direct-capture Perturb-seq Dixit 2016 vs CROP-seq 3'UTR-barcoded Datlinger 2017 vs ECCITE-seq vs Multiome), MOI for sgRNA assignment, escaper-cell filtering (Mixscape, Papalexi 2021), SCEPTRE NB GLM + permutation for low-MOI (Barry 2024 Genome Biol 25:124), the Pertpy framework, factor decomposition, genome-scale Perturb-seq (Replogle 2022 Cell, 2.5M cells), and per-perturbation single-cell DE. Use when running a single-cell CRISPR screen, choosing direct-capture vs CROP-seq architecture, filtering escaper cells, performing single-cell DE, integrating Perturb-seq with pathway analysis, scaling to GW CRISPRi via Replogle protocol, or analyzing multi-omics screens.
tool_type: python
primary_tool: Pertpy
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples checked on pertpy 1.3.0, scanpy 1.12+, anndata 0.13+, pandas 2.2+, numpy 1.26+, scipy 1.12+, sceptre 0.99.0 (R / katsevich-lab/sceptre GitHub HEAD), muon 0.1.9.

**pertpy >= 1.0 breaking changes from the 0.6.x examples some agents may have seen:** `PyDESeq2.test_contrasts()` takes a numeric contrast vector built via `de.contrast(column, baseline, group_to_compare)`, not a `contrast=(column, group, baseline)` tuple; and result columns are `log_fc` / `p_value` / `adj_p_value`, not `log2FoldChange` / `padj`. All code blocks below already use the current API.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show pertpy scanpy anndata`
- R: `packageVersion('sceptre')`; `?sceptre`; `?Seurat::PrepLDA`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Single-Cell Perturb-Seq Analysis

**"Analyze a single-cell pooled CRISPR perturbation screen"** -> Assign sgRNAs to cells, filter unperturbed escapers, normalize counts, fit per-gene differential expression conditioned on perturbation, and rank perturbations by their molecular effect.

- Python: `pertpy` unified framework for Mixscape + SCEPTRE-via-R + differential expression
- R: `sceptre` for low-MOI NB GLM + permutation testing
- Python/R: `Seurat::MixscapeLDA` and downstream

## Experimental Architecture Comparison

| Method | Year | Architecture | Readout | MOI | Single-cell sgRNA detection |
|--------|------|--------------|---------|-----|------------------------------|
| Perturb-seq (Dixit 2016, *Cell*) | 2016 | sgRNA expressed in cassette; direct PCR capture | scRNA-seq | Low-to-moderate (MOI ~0.35-1.4; a minority of cells receive multiple guides, enabling epistasis analysis) | Yes via amplicon-PCR pre-sequencing |
| CROP-seq (Datlinger 2017, *Nat Methods*) | 2017 | hU6-sgRNA cassette placed in the 3' LTR of lentiGuide-Puro; LTR duplication puts the sgRNA in the 3'UTR of the Pol II puromycin-resistance transcript | scRNA-seq | Low (1-2 sgRNAs/cell) | Native via 10X 3' chemistry |
| Perturb-CITE-seq (Frangieh 2021, *Nat Genet*) | 2021 | Adds surface-protein hashtag oligos to CROP-seq | scRNA-seq + ADT (protein) | Low | CROP-seq architecture |
| ECCITE-seq (Mimitou 2019, *Nat Methods*) | 2019 | Surface-protein hashtag with sgRNA-marked cells | scRNA-seq + ADT | Low | Hash + sgRNA |
| Perturb-ATAC (Rubin 2019, *Cell*) | 2019 | scATAC-seq readout | scATAC | Low | sgRNA capture via separate library prep |
| Perturb-multiome (10X) | 2021+ | scRNA + scATAC simultaneously | scRNA + ATAC | Low | Direct capture from sgRNA cassette |
| Replogle GW Perturb-seq (2022, *Cell*) | 2022 | Multiplexed CRISPRi with sgRNA barcoding | scRNA-seq | 1 sgRNA/cell | Direct capture |

**Decision rule:** Standard scRNA + sgRNA at low cost -> CROP-seq. Genome-wide CRISPRi screens -> Replogle's CRISPRi + 10X 3' direct-capture protocol, the gold standard (>2.5M cells; the genome-scale K562 screen targeted ~9,866 expressed genes in Replogle 2022). Protein readout -> Perturb-CITE-seq. Chromatin readout -> Perturb-multiome. Hashed cells + sgRNA -> ECCITE-seq. Low-throughput pilot -> original Dixit Perturb-seq. Genome-scale design and budget: `references/genome-wide-perturb-seq.md`; chromatin readout analysis: `references/multiomic-perturb-seq.md`.

## MOI and sgRNA Assignment

**The central technical challenge:** Each cell must receive exactly one sgRNA (otherwise the perturbation is undefined). At MOI 0.3 (the standard target), ~26% of cells get ≥1 sgRNA, but 4% get ≥2; at MOI 0.5, ~9% of cells get multiple sgRNAs. The cells with multiple sgRNAs must be filtered or analyzed as combinatorial perturbations (see crispr-screens/combinatorial-screens for intentionally-high-MOI paired-guide design).

**Assignment workflow:**

1. **Detect sgRNA reads per cell:** From the sgRNA library prep (direct capture or 3'UTR barcode), count reads per sgRNA per cell.
2. **Threshold:** Most pipelines use 10+ reads of one sgRNA to assign that perturbation.
3. **Multiplets:** Cells with 2+ sgRNAs at >10 reads each are either multi-perturbed (analyzable as combinatorial) or doublets.
4. **Doublet detection:** Use scDblFinder, Scrublet, or AMULET (multiome) to identify doublets independently from sgRNA assignment.

**Goal:** Assign a single perturbation identity (or 'multiplet'/'none') to every cell from the sgRNA counts matrix.

**Approach:** Threshold per-cell sgRNA reads at ≥10 (Pertpy convention); cells exceeding the threshold for exactly one sgRNA are assigned that perturbation; cells with multiple sgRNAs above threshold are flagged as multiplets for filtering or combinatorial analysis.

```bash
python scripts/assign_sgrna.py sgrnas.h5ad --layer sgrna_counts --threshold 10 --out assigned.h5ad
# writes .obs['sgrna_assignment'] = sgRNA name | 'multiplet' | 'none'; or import assign_sgrna() from the script
```

## Escaper Cell Filtering (Mixscape)

**Why this matters:** Not all sgRNA-positive cells actually edit. The escaper fraction is guide- and gene-dependent: Papalexi 2021 measured ~25% escapers for IFNGR2, perturbation rates of 39-92% across four IRF1 guides (i.e. 8-61% escapers), and no detectable perturbation at all for 15 genes. Including escapers dilutes the perturbation effect; Mixscape identifies and filters them.

**Mixscape algorithm:** For each perturbed cell, compute a "perturbation signature" = (its expression) - (mean of K nearest non-targeting-control cells). This signature isolates the perturbation effect from cell-state variation. Cells with perturbation signature similar to NTC distribution are escapers.

```bash
# input: normalized, log1p'd AnnData; NTC label must match your data's actual control label
python scripts/mixscape_filter.py normalized.h5ad --pert-key sgrna_assignment --control NTC --out ko_cells.h5ad
# writes .layers['X_pert'], .obs['mixscape_class'] ('<gene> KO' / '<gene> NP' / control) and
# .obs['mixscape_class_global'] (KO / NP / control); --out holds the KO cells. Seeded (random_state=0)
# so X_pert is reproducible; pertpy >= 1.0.
```

**Critical:** Mixscape can fail when the perturbation has weak phenotype; empirically Mixscape detects perturbations with log-fold-change <-0.5 (depletion) reliably, but weaker effects collapse into the NTC distribution. For genome-wide screens, run Mixscape per perturbation; for low-effect perturbations, trust the assignment without filtering.

## SCEPTRE for Low-MOI Differential Expression

SCEPTRE (NB GLM + conditional resampling, calibrated FDR) is the low-MOI DE method of choice, run in R. Full method and code: `references/sceptre-low-moi.md`.

## Pertpy Unified Framework

**Pertpy** (https://pertpy.readthedocs.io) integrates Mixscape, distance-based perturbation comparison, EdgeR/PyDESeq2/WilcoxonTest DE, and factor models in a single AnnData-based interface. For SCEPTRE specifically, invoke the R sceptre package separately (Pertpy does not wrap it).

The end-to-end recipe (load papalexi_2021, merge `gene_target` from the MuData, keep raw counts in `layers['counts']`,
Mixscape with `control='NT'`, per-perturbation PyDESeq2 contrasts) is `examples/run_pertpy.py`. It defaults to a bounded
600-cell / 2,000-gene smoke run across two perturbations; provide `--full` only with a separately planned tens-of-GiB memory
budget. The DE call it makes:

```python
de = pt.tl.PyDESeq2(adata_ko, design='~gene_target', layer='counts')   # raw counts, not log-normalized
de.fit()
results_df = de.test_contrasts(de.contrast('gene_target', 'NT', 'GENE_X'))   # log_fc / p_value / adj_p_value
```

## Genome-Wide Perturb-Seq (Replogle 2022)

Genome-scale CRISPRi design, cell/channel budget and the scaling formula for a different gene count: `references/genome-wide-perturb-seq.md`.

## Factor-Based Analysis

Shared-factor decomposition of perturbation effects (FR-Perturb, standalone CLI): `references/factor-decomposition.md`.

## Multiomic Perturb-seq (RNA + ATAC)

RNA + ATAC Perturb-seq (10X Multiome): propagate the Mixscape call to the ATAC modality and test differential accessibility: `references/multiomic-perturb-seq.md`.

## Reference Files

| File | Read when |
|------|-----------|
| `references/sceptre-low-moi.md` | Calibrated single-cell DE for a low-MOI screen (R `sceptre` pipeline) |
| `references/genome-wide-perturb-seq.md` | Designing or budgeting a genome-scale CRISPRi Perturb-seq (Replogle 2022) |
| `references/factor-decomposition.md` | Decomposing perturbation effects into shared factors (FR-Perturb) |
| `references/multiomic-perturb-seq.md` | RNA + ATAC (10X Multiome) Perturb-seq: differential accessibility per perturbation |
| `references/optional-methods.md` | Installing or falling back from optional FR-Perturb and Seurat methods |

## Failure Modes

### Low sgRNA detection per cell

**Trigger:** Direct-capture method on CROP-seq library, or 3'UTR barcoding on direct-capture library.
**Mechanism:** Architecture mismatch -- the sgRNA can't be detected by the wrong library prep.
**Symptom:** sgRNA assignment rate <50% of cells.
**Fix:** Match library prep to architecture; for CROP-seq, use 10X 3' chemistry; for direct-capture Perturb-seq, use the Dixit amplicon-PCR pre-sequencing.

### Mixscape filters too many cells as escapers

**Trigger:** Weak perturbation phenotype; Mixscape's NTC-subtracted signature is similar to NTC null.
**Mechanism:** Mixscape assumes a detectable signal; weak knockdown is misclassified as escaper.
**Symptom:** >50% of perturbed cells classified as "NP" (non-perturbed); known essentials show no effect.
**Fix:** Lower Mixscape stringency; skip Mixscape for low-effect perturbations; verify Cas9 expression first.

### Doublet contamination drives apparent multi-perturbation cells

**Trigger:** High cell density loading on 10X channels.
**Mechanism:** Two cells in one droplet appear to carry two sgRNAs.
**Symptom:** "Multiplet" rate >5% after sgRNA assignment.
**Fix:** Reduce cell loading per channel (5,000-7,000 instead of 10,000); run Scrublet or scDblFinder; remove doublets before sgRNA assignment.

### MAST or Wilcoxon over-call hits

**Trigger:** Using parametric DE tools on sparse, zero-inflated scRNA-seq.
**Mechanism:** These tools assume Gaussian or simpler null; single-cell data has zero-inflation that makes them over-confident.
**Symptom:** Thousands of significant DE genes per perturbation; FDR uncalibrated.
**Fix:** Use SCEPTRE (permutation-based NB GLM); Barry 2024 benchmark shows this is the only method with calibrated FDR.

### Genome-scale Perturb-seq with insufficient cells per perturbation

**Trigger:** <500 cells per perturbation in genome-scale experiment.
**Mechanism:** DE estimation requires sufficient cells per condition; <500 lacks power for moderate effects.
**Symptom:** Inconsistent hit calls across replicates; pathway analysis non-specific.
**Fix:** Scale up cell numbers; or run focused (sub-genome) Perturb-seq with more cells per pert.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| MOI for single sgRNA per cell | 0.3 | Poisson math; ~26% infected, 4% multi-infected |
| sgRNA assignment threshold | ≥10 reads of one sgRNA | Pertpy / direct-capture convention |
| Multiplet rate (post-doublet filter) | <5% | Typical 10X 3' chemistry |
| Mixscape KO retention | Guide-dependent; 39-92% observed | Papalexi 2021 |
| Cells per perturbation (DE power) | 500-1,000 minimum genome-scale; 1,000-2,000 focused (specific module); 5,000+ single-pert deep; 2,000+ per pair combinatorial | Power convention (Replogle 2022 screened at a median >100) |
| SCEPTRE permutations | 1,000+ | Barry 2024 |
| Genes per cell (QC) | ≥500-1,000 | Standard scRNA QC |
| Mt% threshold | <15-20% | Standard scRNA QC |
| Doublet detection threshold | scDblFinder, Scrublet defaults | Methods agree |
| NTC (non-targeting control) representation | ~5% of library | Standard pooled-screen library design |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Low sgRNA detection | Architecture mismatch | Match library prep |
| Too many escapers in Mixscape | Weak phenotype | Skip Mixscape; verify Cas9 |
| Inflated DE hits | MAST / Wilcoxon used | Switch to SCEPTRE |
| Inconsistent gene effects between channels | Channel batch effect | Add channel as covariate in SCEPTRE |
| Multiplet rate >10% | Over-loading cells | Reduce loading; doublet filter |
| Per-pert DE with <100 cells | Insufficient power | Increase cell numbers; or accept low resolution |

## Scope

This Skill analyzes single-cell pooled CRISPR screen data for research purposes. A perturbation's molecular effect here (a DE gene, an escaper-filtered phenotype, a factor loading) is a research finding, not a validated therapeutic target or a patient-treatment recommendation -- a screen-level knockdown effect does not by itself establish clinical benefit or harm for a patient. Do not use this Skill's output to make or imply a diagnostic or treatment decision for an individual.

## References

- Dixit A et al. 2016. *Cell* 167:1853. Original Perturb-seq.
- Datlinger P et al. 2017. *Nat Methods* 14:297. CROP-seq.
- Frangieh CJ et al. 2021. *Nat Genet* 53:332. Perturb-CITE-seq.
- Mimitou EP et al. 2019. *Nat Methods* 16:409. ECCITE-seq.
- Rubin AJ et al. 2019. *Cell* 176:361. Perturb-ATAC.
- Papalexi E et al. 2021. *Nat Genet* 53:322. Mixscape.
- Barry T, Mason K, Roeder K, Katsevich E. 2024. *Genome Biol* 25:124. SCEPTRE for low-MOI Perturb-seq.
- Replogle JM et al. 2022. *Cell* 185:2559. Genome-wide Perturb-seq.
- Heumos L et al. 2026. *Nat Methods* 23:350-359. DOI 10.1038/s41592-025-02909-7. Pertpy framework.
- Jiang L et al. 2025. *Nat Cell Biol* 27:505. Mixscale (perturbation-strength-aware Perturb-seq).

## Related Skills

- crispr-screens/library-design - Direct-capture vs CROP-seq library design
- crispr-screens/screen-qc - sgRNA assignment rates as QC
- crispr-screens/mageck-analysis - Pseudobulk analysis as alternative
- crispr-screens/hit-calling - Pseudo-bulk hit calling alternative
- single-cell/preprocessing - scRNA-seq preprocessing
- single-cell/clustering - Post-DE clustering
- single-cell/multimodal-integration - Multiome Perturb-seq
- single-cell/perturb-seq - General single-cell screen analysis
- pathway-analysis/go-enrichment - Pathway enrichment of perturbation hits
