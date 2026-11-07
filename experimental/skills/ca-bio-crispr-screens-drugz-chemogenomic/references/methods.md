## Why drugZ for Drug Screens (not MAGeCK)

| Property | drugZ | MAGeCK RRA | MAGeCK MLE |
|----------|-------|------------|-------------|
| Bidirectional sensitivity | YES (sensitizer + resistance same scale) | Asymmetric (neg/pos separately) | Asymmetric |
| Drug-anchored baseline | YES (drug vs vehicle) | Either (drug vs vehicle or vs Day 0) | Either |
| Sensitivity to small effects | Highest (bidirectional Z; Colic et al. 2019) | Moderate | Moderate |
| Statistical framework | Empirical-Bayes windowed Z-score on guide-level log fold change | NB + alpha-RRA | NB GLM with design matrix |
| Handles guide-level noise | sgRNA-level z aggregation | Rank-based aggregation | Built-in guide-efficacy term (optional) |
| Best for | Drug-modifier / chemogenomic screens | General essentiality / standard 2-condition | Time course / multi-condition |

**Why MAGeCK is suboptimal for drug screens:** MAGeCK's RRA was designed for two-condition essentiality; drug-vs-vehicle screens often have small effect sizes (10-30% sgRNA shift) that RRA rank-based aggregation under-detects. drugZ uses parametric Z-scoring tuned for these small effects.

**Benchmark (Colic et al. 2019):** On DNA-damage-response chemogenomic screens, drugZ hits were far more strongly enriched for the expected pathway (DDR) than STARS, MAGeCK, edgeR or RIGER hits across FDR thresholds, reflecting better sensitivity to the moderate fitness defects typical of drug-gene interactions. Compare methods on expected-pathway enrichment, not raw hit count.

## The drugZ Algorithm (under the hood)

1. For each sgRNA, compute log2-fold-change drug vs vehicle: `LFC_drug_vs_veh`
2. Compute an empirical-Bayes Z per sgRNA: `Z = LFC / eb_std`, where `eb_std` is the standard deviation of a sliding window of guides with similar control abundance (`--half_window_size`, default 500), smoothed monotonically
3. Per gene, sum Z across all sgRNAs targeting it: `sumZ = sum(Z_sgRNA)`
4. Normalize and re-standardize across genes: `normZ = zscore(sumZ / sqrt(numObs))`
5. Compute a one-sided p-value per direction: synth (sensitizer = negative normZ) and supp (resistance = positive normZ)
6. Benjamini-Hochberg FDR correction per direction

**Critical:** Vehicle vs drug, NOT Day 0 vs drug. Day-0 baseline conflates proliferation effects with drug effects.

**Library size vs `--half_window_size`:** the sliding window in step 2 is indexed absolutely, so a
count table with fewer than about 4x `--half_window_size` guides raises
`IndexError: single positional indexer is out-of-bounds` instead of a useful message. The default 500
therefore needs a genome-scale library (tens of thousands of guides). For a small pilot or arrayed
sub-library, set `--half_window_size` to roughly a quarter of the total guide count, and treat the
per-guide variance estimate as approximate. This is a *total guide count* constraint, separate from
the 4-6 sgRNAs *per gene* needed for a stable per-gene Z.


## Vehicle vs Day-0 Reference: Critical Decision

**Why this matters:** Drug screen analysis can compare drug to:
1. **Vehicle (DMSO / carrier)** -- isolates drug-specific effect; correct anchor.
2. **Day 0 (initial library)** -- conflates proliferation, drug, and vehicle effects.

```
counts at Day 0          (no perturbation; cloning baseline)
    |
    v
counts at Day 7 - Vehicle (proliferation only; what survives in normal culture)
counts at Day 7 - Drug    (proliferation + drug effect)
    |
    v
Drug effect = LFC(Drug vs Vehicle)         # CORRECT
Wrong:       LFC(Drug vs Day 0)            # confounds drug with general proliferation
```

drugZ specifically requires `-c` to name the vehicle samples. Always include matched vehicle controls in drug screens.

## Drug-Dose and Time-Course Designs

**drugZ for dose-response:** Not natively designed for dose; instead, run drugZ separately at each dose vs vehicle, then look for genes with consistent direction across doses.

```bash
for DOSE in low mid high; do
    python drugz.py \
        -i counts.txt \
        -o drugz_${DOSE}.txt \
        -c Veh_r1,Veh_r2 \
        -x Drug${DOSE}_r1,Drug${DOSE}_r2
done
```

Then aggregate. **Dose consistency rule:** a dose-consistent hit keeps the same `normZ` sign at every
tested dose and reaches FDR < 0.05 at the highest dose; report the others as dose-inconsistent rather
than dropping them silently. `|normZ|` growing with dose is supporting evidence, not a requirement
(saturation at the top dose is common); `scripts/dose_consistent_hits.py` applies the rule:

```bash
# DOSE=FILE pairs, low to high; --top-dose names the highest dose. Also importable: from dose_consistent_hits import dose_consistent_hits
python scripts/dose_consistent_hits.py --top-dose high low=drugz_low.txt mid=drugz_mid.txt high=drugz_high.txt
# add --direction supp for suppressors, --fdr 0.05 (default), --out hits.tsv
```

**For multi-condition drug-screens** (time × drug × cell-line), use MAGeCK MLE with explicit design matrix instead -- MLE handles multi-factorial; drugZ does not.

## Comparison: drugZ vs MAGeCK MLE for Drug Screen

**Goal:** When to use each method.

| Question | drugZ | MAGeCK MLE |
|----------|-------|-------------|
| Single drug, single dose, vehicle vs drug | YES (preferred) | Acceptable |
| Multiple doses, drug response curve | Per-dose drugZ + meta | YES (preferred with dose covariate) |
| Time course at single dose | Per-timepoint drugZ + meta | YES (preferred with time covariate) |
| Drug + cell-line panel | Per-line drugZ + meta | YES (or Chronos) |
| Combinatorial drug pairs | Per-pair drugZ + meta | YES (preferred with interaction) |
| Synergy / antagonism detection | Limited (per-drug calling only) | YES (interaction term in MLE) |
| Small effect sizes (LFC <0.5) | Highest sensitivity | Lower sensitivity |
| Heavy selection (>40% guides change) | OK | Norm needs control sgRNAs |
| Essentiality plus drug effect | drugZ for the drug effect | BAGEL2 for essentiality |

Dose consistency: see the rule in "Drug-Dose and Time-Course Designs". Failure symptoms (no hits, essentials dominating, unstable or dose-blind hits, crashes): `references/failure-modes.md`.

**Reconciliation:** For simple drug-modifier screens with one drug and one vehicle, run both drugZ and MAGeCK MLE; hits called by both are high confidence; drugZ-only hits at low LFC need orthogonal validation (drug + arrayed validation).


## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Sensitizer hit | `fdr_synth < 0.05` | Colic et al. 2019; BH-corrected |
| Suppressor hit | `fdr_supp < 0.05` | Same |
| High-confidence sensitizer | `fdr_synth < 0.01 AND normZ < -3` | Conservative |
| Pseudocount default | 5 | Colic et al. 2019 |
| Min sgRNAs per gene for stable Z | 4-6 | Below this, Z varies between runs |
| Vehicle replicates needed | 3+ | For stable Z null distribution |
| Drug replicates needed | 3+ | For per-gene sumZ stability |
| Replicate Pearson within an arm | > 0.85 | Pre-drugZ QC; below it, see "Low replicate concordance" in `references/failure-modes.md` |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| No hits | Wrong control samples (Day-0 instead of vehicle) | Re-run with vehicle |
| Hits dominated by essentials | Essentials inflate null | Use `-r` with a comma-list of CEGv2 (`references/reference-gene-removal.md`) |
| Unstable hits across runs | Too few sgRNAs/gene | Use 6+ sgRNAs/gene library |
| Drug-target appears in suppressor | Real biology | Annotate separately |
| MAGeCK and drugZ disagree | Different statistical sensitivity | drugZ more sensitive; trust for chemogenomic |
| Inconsistent between doses | Real dose effect | Require consistency across doses |
| `IndexError: single positional indexer is out-of-bounds` | `--half_window_size` too large for the total guide count (needs about 4x) | Set `--half_window_size` to about 1/4 of total guides; see "Library size vs `--half_window_size`" |
