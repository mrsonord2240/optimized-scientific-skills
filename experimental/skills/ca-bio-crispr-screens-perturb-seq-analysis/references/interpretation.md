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
