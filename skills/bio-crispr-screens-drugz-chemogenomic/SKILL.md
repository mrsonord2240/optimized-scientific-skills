---
name: bio-crispr-screens-drugz-chemogenomic
description: Analyzes CRISPR drug-modifier (chemogenomic) screens with drugZ (Colic et al. 2019 Genome Med), a bidirectional Z-score method that identifies synthetic-lethal sensitizing genes and resistance-conferring suppressor genes from vehicle vs drug comparisons. Covers vehicle-anchored design (not Day-0), the bidirectional Z math giving greater sensitivity to small-effect hits than MAGeCK / STARS / edgeR / RIGER on drug screens, per-gene sumZ and normZ, synth (sensitizer) vs supp (suppressor) FDR, multi-dose handling, integration with control sgRNAs, and comparison with MAGeCK MLE with dose covariate. Use when running a drug-modifier CRISPR screen, identifying sensitizing or resistance genes for a drug candidate, choosing drugZ vs MAGeCK MLE for chemogenomic analysis, troubleshooting low-effect drug screens where MAGeCK lacks sensitivity, or designing a drug-screen layout (vehicle vs drug arms).
tool_type: cli
primary_tool: drugZ
license: MIT
---

## Version Compatibility

Reference examples tested with: drugZ Aug-2019+ (hart-lab/drugz; checked 2026-09-16 against the current master `drugz.py`), MAGeCK 0.5.9+, pandas 2.2+, numpy 1.26+, scipy 1.12+, statsmodels 0.14+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `python drugz.py --help` (the repo has no setup.py, so there is no `drugz` console script)
- GitHub: install via `git clone https://github.com/hart-lab/drugz`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## drugZ Chemogenomic Analysis

**"Identify genes that sensitize or confer resistance to my drug in a CRISPR screen"** -> Compare drug-treated vs vehicle-treated arms (NOT Day-0 baseline) using bidirectional Z-scores per sgRNA, sum to per-gene normalized Z, and rank genes for sensitizer (synthetic lethal) vs suppressor (resistance) phenotype.

- CLI: `python drugz.py -i counts.txt -o drugz.txt -c Vehicle_r1,Vehicle_r2 -x Drug_r1,Drug_r2`
- Python: programmatic via `drugz.drugZ_analysis(args)` (takes an argparse Namespace)
- Workflow: vehicle-anchored counts -> Z-scoring -> per-gene summation -> direction-specific FDR

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

## Run drugZ on a Drug-Modifier Screen

**Goal:** Quantify per-gene sensitizing and suppressor effects from a chemogenomic screen.

**Approach:** Run `drugz.py` with vehicle and drug sample columns; output per-gene sumZ, normZ, and direction-specific p-values + FDR.

```bash
git clone https://github.com/hart-lab/drugz
cd drugz

# Standard drug screen comparison:
# Vehicle (DMSO or carrier) replicates: Veh_r1, Veh_r2, Veh_r3
# Drug-treated replicates: Drug_r1, Drug_r2, Drug_r3

python drugz.py \
    -i counts.txt \                       # input read-count file (tab-separated)
    -o drugz_output.txt \                  # output file
    -c Veh_r1,Veh_r2,Veh_r3 \              # control samples (comma-separated)
    -x Drug_r1,Drug_r2,Drug_r3 \           # treated samples (comma-separated)
    -r RPS3,RPL11,EIF3A \                  # OPTIONAL: comma-delimited GENE NAMES to exclude (not a file path)
    -p 5                                   # pseudocount (default 5)

# Output: drugz_output.txt with columns:
#   GENE, sumZ, numObs, normZ, pval_synth, rank_synth, fdr_synth, pval_supp, rank_supp, fdr_supp
```

**Output columns:**

| Column | Meaning |
|--------|---------|
| `GENE` | Gene symbol |
| `numObs` | Number of non-zero guide x replicate observations |
| `sumZ` | Summed per-sgRNA Z-score |
| `normZ` | sumZ / sqrt(numObs), re-standardized across genes |
| `pval_synth` | One-sided p-value for sensitizer (negative effect; gene KO sensitizes to drug) |
| `rank_synth` | Rank for sensitizers |
| `fdr_synth` | BH-corrected FDR for sensitizers |
| `pval_supp` | One-sided p-value for suppressor (positive effect; gene KO confers resistance) |
| `rank_supp` | Rank for suppressors |
| `fdr_supp` | BH-corrected FDR for suppressors |

**Interpretation:**
- Sensitizers (synthetic lethal): `fdr_synth < 0.05` -- loss of these genes makes cells more sensitive to drug. Examples: PARPi targets BRCA1/2; cisplatin sensitizes ERCC.
- Suppressors (resistance): `fdr_supp < 0.05` -- loss of these genes confers resistance. Examples: drug-efflux genes; drug target itself paradoxically.

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

Then aggregate. A dose-consistent hit is one that keeps the same sign at every tested dose and
reaches FDR < 0.05 at the highest dose; report the others as dose-inconsistent rather than dropping
them silently:

```python
import pandas as pd

def dose_consistent_hits(dose_files, top_dose, fdr=0.05, direction='synth'):
    """dose_files: {'low': 'drugz_low.txt', 'mid': ..., 'high': ...}; top_dose: key of the highest dose.

    A hit is dose-consistent if normZ has the same sign at every dose and passes FDR at the top dose.
    """
    frames = {d: pd.read_csv(f, sep='\t').set_index('GENE') for d, f in dose_files.items()}
    normz = pd.DataFrame({d: f['normZ'] for d, f in frames.items()}).dropna()
    sign_ok = (normz.gt(0).all(axis=1)) | (normz.lt(0).all(axis=1))
    top = frames[top_dose]
    passes = top['fdr_%s' % direction] < fdr
    hits = normz[sign_ok & passes.reindex(normz.index).fillna(False)].copy()
    hits['normZ_top_dose'] = top['normZ'].reindex(hits.index)
    hits['monotonic'] = (normz.abs().diff(axis=1).iloc[:, 1:] >= 0).all(axis=1)  # |normZ| grows with dose
    return hits.sort_values('normZ_top_dose')
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

**Dose consistency rule:** same sign of `normZ` at every tested dose and FDR < 0.05 at the highest dose. Report `|normZ|` increasing with dose as supporting evidence, not as a requirement -- saturation at the top dose is common.

**Reconciliation:** For simple drug-modifier screens with one drug and one vehicle, run both drugZ and MAGeCK MLE; hits called by both are high confidence; drugZ-only hits at low LFC need orthogonal validation (drug + arrayed validation).

## Removing Reference Genes from the Analysis

**Goal:** Keep reference essential or control genes from inflating the Z-score null distribution.

**Approach:** `-r` takes a **comma-delimited list of gene symbols, not a file path** (`drugz.py` does
`args.remove_genes.split(',')`). Passing a filename removes nothing and still exits 0. The named genes
are dropped from the input before any Z-scoring, so they also disappear from the output file: `-r` is
an exclusion, not a re-weighting.

```bash
# Inline list
python drugz.py \
    -i counts.txt \
    -o drugz_clean.txt \
    -c Veh_r1,Veh_r2 \
    -x Drug_r1,Drug_r2 \
    -r RPS3,RPL11,EIF3A,POLR2A,CDK1
```

```bash
# From a reference set: take the first column, skip the header, join with commas.
# CEGv2.txt (hart-lab/bagel) is tab-separated with a header (GENE, HGNC_ID, ENTREZ_ID);
# joining its raw lines gives tokens like "AARS<TAB>HGNC:20<TAB>16", which match no gene
# and silently exclude nothing.
CEG=$(tail -n +2 CEGv2.txt | cut -f1 | paste -sd, -)
python drugz.py -i counts.txt -o drugz_clean.txt -c Veh_r1,Veh_r2 -x Drug_r1,Drug_r2 -r "$CEG"
```

**Verify the exclusion happened** -- the tool cannot tell you it matched nothing:

```bash
# genes in the unfiltered output but not the filtered one; should equal the number you excluded
comm -23 <(cut -f1 drugz_output.txt | sort) <(cut -f1 drugz_clean.txt | sort) | wc -l
```

**When to use:** If pilot drugZ runs show many essential genes appearing as "sensitizers" purely
because they drop out under any condition, removing them gives a cleaner drug-specific signal.

## Failure Modes

### drugZ shows no synthetic-lethal hits despite known sensitizing genes

**Trigger:** Comparing drug vs Day-0 instead of drug vs vehicle.
**Mechanism:** Day-0 comparison conflates drug effect with normal-culture proliferation; essential genes drop in both conditions, masking drug-specific sensitization.
**Symptom:** PARPi screen shows no sensitization at BRCA1/BRCA2 despite expected biology.
**Fix:** Re-run with vehicle samples passed to `-c`. The drug-vs-vehicle is the canonical comparison.

### High false-positive rate among essential genes

**Trigger:** Essential genes drop out in both vehicle and drug arms; small relative shift gives misleadingly high Z.
**Mechanism:** drugZ's Z-score is symmetric; essential genes drop in both arms but slightly more in drug -> "synthetic lethal" call.
**Symptom:** Hit list dominated by RPS, RPL, EIF essentials.
**Fix:** Use `-r` with a comma-delimited list of essential gene names to exclude (see "Removing Reference Genes"), then check that the excluded genes really are missing from the output; or filter the output post-hoc.

### Unstable hits across libraries or sub-samples

**Trigger:** Insufficient sgRNAs per gene; small effect sizes.
**Mechanism:** drugZ's per-gene sumZ depends on enough sgRNAs to be stable; with 3-4 sgRNAs/gene, single-guide noise drives the ranking.
**Symptom:** Top hits move when you re-sequence, sub-sample replicates, or switch library.
**Note:** Re-running drugZ on the same input cannot show this. drugZ has no sampling step and no seed: identical input gives a byte-identical output file, so a rerun is not a stability check.
**Fix:** Use a 6+ sgRNAs/gene library (Avana, Dolcetto); check stability by holding out a replicate or bootstrapping the guides yourself; or use MAGeCK MLE.

### drugZ ignores dose information

**Trigger:** Multi-dose screen analyzed at highest dose only.
**Mechanism:** drugZ doesn't model dose; running at one dose loses the dose-response information.
**Symptom:** Hits at high dose may be dose-specific (not true responders).
**Fix:** Run drugZ at each dose; require consistency across doses for high-confidence hits.

### Drug-target gene appears as "suppressor"

**Trigger:** Loss of drug target reduces drug binding, increasing drug resistance.
**Mechanism:** Real biology -- drug target itself is a resistance gene from a KO perspective.
**Symptom:** Drug-target gene like PARP1 appears in suppressor list for PARPi screen.
**Fix:** Expected biology. Annotate the drug target separately. The suppressor list is correct.

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

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| No hits | Wrong control samples (Day-0 instead of vehicle) | Re-run with vehicle |
| Hits dominated by essentials | Essentials inflate null | Use `-r` with a comma-list of CEGv2 |
| Unstable hits across runs | Too few sgRNAs/gene | Use 6+ sgRNAs/gene library |
| Drug-target appears in suppressor | Real biology | Annotate separately |
| MAGeCK and drugZ disagree | Different statistical sensitivity | drugZ more sensitive; trust for chemogenomic |
| Inconsistent between doses | Real dose effect | Require consistency across doses |

## References

- Colic M et al. 2019. *Genome Medicine* 11:52. drugZ algorithm and chemogenomic-interaction benchmark.
- Olivieri M et al. 2020. *Cell* 182:481. DDR chemogenomic screens with drugZ.
- Behan FM et al. 2019. *Nature* 568:511. Project Score; genome-wide cancer-dependency screens for target prioritization.

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK MLE alternative for multi-condition drug screens
- crispr-screens/bagel-essentiality - BAGEL2 alternative; sensitive to tumor-suppressor / drug-target
- crispr-screens/hit-calling - Cross-method decision tree including drugZ
- crispr-screens/screen-qc - Pre-drugZ QC including replicate concordance
- crispr-screens/library-design - 6+ sgRNAs/gene library for stable Z
- crispr-screens/copy-number-correction - Pre-correction for cancer-line drug screens
- crispr-screens/base-editing-analysis - Variant-function drug-modifier screens
- pathway-analysis/go-enrichment - Functional analysis of drug-modifier hits
- clinical-databases/clinvar-lookup - Clinical interpretation of drug targets
