---
name: bio-crispr-screens-copy-number-correction
category: Data Analysis
description: Corrects the gene-independent copy-number artifact in CRISPR-Cas9 screens (Aguirre 2016 / Munoz 2016 Cancer Discov) where amplified loci appear essential from DNA-damage burden of simultaneous cuts. Covers the gene-independent DNA-damage / G2-arrest mechanism, CRISPRcleanR (Iorio 2018) unsupervised pre-hoc correction, CERES (Meyers 2017) joint CN + gene-effect model, Chronos (Dempster 2021) DepMap-standard population-dynamics + CN model with lowest residual bias, the decision tree by data availability, the Spearman LFC-vs-CN diagnostic, focal-amplification examples (ERBB2 in HER2+, MYC in colorectal, FGFR1 in head and neck), and CRISPRi/a alternatives that bypass the artifact. Use when screening cancer cell lines, diagnosing essentiality at amplified loci, choosing CRISPRcleanR / CERES / Chronos, deciding whether CN correction is needed before MAGeCK / BAGEL2 / drugZ, or switching from Cas9 to CRISPRi.
tool_type: mixed
primary_tool: CRISPRcleanR
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: CRISPRcleanR 3.0+ (R; github.com/francescojm/CRISPRcleanR), Chronos 2.3.15 (`crispr_chronos`; checked 2026-09-16), CERES (legacy, superseded by Chronos), pandas 2.2+, numpy 1.26+, scipy 1.12+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('CRISPRcleanR')`; `?ccr.GWclean`
- Python: `pip show crispr_chronos`; `python -c 'import chronos; print(chronos.__file__)'`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Copy-Number Bias Correction in CRISPR Screens

**Scope:** this Skill analyses screen data. It does not support treatment, therapy or other clinical
recommendations for an individual; an amplification that survives correction is a research finding,
not evidence for a therapy choice. Decline the clinical part of such a request and answer the
computational part.

**"Correct copy-number artifacts in my cancer-cell-line screen"** -> Identify gene-independent depletion at amplified loci, apply CRISPRcleanR (pre-hoc, unsupervised, position-based) or Chronos (joint model, supervised with CN profile) to remove the artifact, then proceed to hit calling on corrected data.

- R: `CRISPRcleanR::ccr.GWclean()` for unsupervised pre-hoc correction (no CN profile required)
- Python: Chronos (`crispr_chronos`) for joint cell-population dynamics + CN modeling
- Python: CERES (legacy, superseded by Chronos)

## The Copy-Number Artifact (Mechanism)

**Aguirre AJ et al 2016 *Cancer Discov* 6:914** and **Munoz DM et al 2016 *Cancer Discov* 6:900** demonstrated that focal amplification regions in cancer cell lines appear systematically "essential" in CRISPR-Cas9 screens, independent of the gene's actual biology. The mechanism:

1. A focal amplification creates 4-50+ copies of a genomic region.
2. Each sgRNA targeting a gene in that region cuts at all copies simultaneously.
3. Multiple cuts trigger a DNA-damage response and G2 arrest, in both TP53-mutant and TP53-wild-type lines but with larger magnitude in wild-type (Aguirre 2016).
4. Cells arrest in G2 phase; the sgRNA appears depleted because its bearer cells don't proliferate.
5. The depletion is proportional to the number of simultaneous cuts, not the gene's essentiality.

**Consequence:** ERBB2 appears essential in HER2-amplified SK-BR-3. MYC appears essential in MYC-amplified colorectal lines (10+ copies). FGFR1 appears essential in FGFR1-amplified head-and-neck lines. These are all false positives.

**Affects:** All Cas9-KO screens in cancer cell lines. Universal, not conditional. Cannot be remediated by sequencing depth, library size, or replicate count. Requires explicit correction.

The p53-dependence of Cas9-cut toxicity in general was characterized later, by Haapaniemi 2018 and Ihry 2018.

**Bypassed by:**
- CRISPRi (catalytically dead Cas9, no DNA damage) -> no artifact
- CRISPRa (catalytically dead Cas9) -> no artifact
- Base editing (single-strand nick + deaminase) -> reduced artifact
- Prime editing (nick + RT) -> reduced artifact

## Correction Method Decision Tree

| Available data | Recommended method | Why |
|----------------|---------------------|-----|
| Cell-line panel without matched CN profile | CRISPRcleanR | Unsupervised; uses genomic position only |
| Single cell line with matched WGS/SNP-array CN | CRISPRcleanR | Chronos' CN correction (`alternate_CN`) refuses fewer than 3 cell lines, matched CN or not; it says so and points at CRISPRcleanR |
| >=3 cell lines, CN available | Chronos | `alternate_CN`'s minimum; it fits the CN-effect curve across lines |
| DepMap-scale (1000+ cell lines, longitudinal) | Chronos | Population-dynamics + screen quality + CN; DepMap quarterly standard |
| Single cell line, multi-timepoint | Chronos for gene effects, CRISPRcleanR for CN correction | Chronos trains and scores fine on one line, but cannot CN-correct below 3 lines |
| Need to integrate with downstream MAGeCK | CRISPRcleanR (pre-hoc) | Outputs corrected counts for any downstream tool |
| Multiple cell lines + multiple batches | Chronos | Joint modeling of all dimensions |

## CRISPRcleanR (Iorio 2018) - Unsupervised Pre-Hoc

**Goal:** Correct copy-number bias without requiring matched CN profile by detecting position-based systematic enrichment / depletion patterns.

**Approach:** Order sgRNAs by chromosomal coordinate; detect segments where sgRNAs show systematic depletion (or enrichment) inconsistent with single-gene biology; shift these segments toward the global mean. The intuition: focal amplifications create depletion bands extending tens to hundreds of kb; non-amplified essential genes are punctate.

```r
library(CRISPRcleanR)

# Load library annotation (sgRNA -> chromosomal coordinates)
data(KY_Library_v1.0)   # KY library; replace with your library annotation
# OR use ccr.PrepareAnnotations() to make custom

# Load count data with first 2 cols: sgRNA, gene, then sample counts
counts <- read.table('counts.txt', header=TRUE, sep='\t')

# 1. Normalize and compute logFC
norm_counts <- ccr.NormfoldChanges(filename='counts.txt', min_reads=30,
                                     EXPname='my_screen',
                                     libraryAnnotation=KY_Library_v1.0)

# 2. Compute genome-sorted sgRNA fold changes
gw_log_fc <- ccr.logFCs2chromPos(norm_counts$logFCs,
                                   KY_Library_v1.0)

# 3. Apply CRISPRcleanR correction
corrected <- ccr.GWclean(gw_log_fc, display=TRUE, label='my_screen')
# Output: corrected$corrected_logFCs and corrected$segments
# corrected_logFCs can replace LFCs downstream

# 4. Re-derive corrected counts for downstream MAGeCK
corrected_counts <- ccr.correctCounts('my_screen',
                                        norm_counts$norm_counts,
                                        corrected,
                                        KY_Library_v1.0,
                                        OutDir='./')
```

**Key parameter:** `min_reads=30` is the lower-count threshold for inclusion. This must match the library-coverage strategy; too high removes legitimate guides, too low keeps noisy guides.

**Output:** Pre-corrected LFCs and counts that can be fed into MAGeCK / BAGEL2 / drugZ as if they were the original screen data. The correction is independent of CN profile (unsupervised) and works on cell lines without matched WGS.

## Chronos (Dempster 2021) - Joint Population-Dynamics + CN Model

**Goal:** Estimate gene fitness while jointly accounting for copy-number-driven depletion, screen quality, and longitudinal cell-population dynamics.

**Approach:** Model the cell population over time as an ODE driven by per-gene fitness effects; add a separate term for copy-number-driven depletion; estimate all parameters via maximum-likelihood with regularization. Outputs a "gene effect score" normalized against the empirical distributions of essential and non-essential reference genes.

```python
# Chronos (pip install crispr_chronos, or pip install git+https://github.com/broadinstitute/chronos)
import numpy as np
import chronos
from chronos.hit_calling import get_probability_dependent

# Chronos initializes model parameters with NumPy draws. Seed *before* constructing the model,
# and record the Chronos version, seed, input hashes, epochs, and platform in the run manifest.
# This makes a CPU rerun materially comparable; do not promise bitwise-identical estimates across
# TensorFlow/BLAS versions or hardware. Re-run important panels and report their stability.
np.random.seed(20260924)  # replace with and record a project-specific seed

# Every input is a dict keyed by library name, not a bare DataFrame, and Chronos checks the
# schema below (chronos.check_inputs) before training:
#
# 1. readcounts    rows = sequence_ID (one per sequenced sample), columns = sgRNA.
#                  This is the transpose of a MAGeCK count table; Chronos rejects the other
#                  orientation with "Chronos expects readcounts to have guides as columns,
#                  sequence IDs as rows. Is your data transposed?"
# 2. sequence_map  columns sequence_ID, cell_line_name, days, pDNA_batch. Plasmid samples get
#                  cell_line_name 'pDNA' and days 0, and every pDNA_batch used by a late
#                  timepoint needs a pDNA row in the same library.
# 3. guide_gene_map  columns sgrna, gene; one gene per sgRNA (duplicated sgRNAs are rejected).
# 4. negative_control_sgrnas  dict library -> sgRNAs of non-targeting or known non-essential
#                  genes. Optional in the signature only: without it, chronos.Chronos(...)
#                  construction itself raises ValueError('excess_variance was passed as dict
#                  without key for <library>: {}') -- before train() is ever called. An
#                  UnboundLocalError for 'prior_variance' does occur one frame deeper, inside
#                  _estimate_excess_variance, but Chronos catches it internally and re-raises
#                  the ValueError above every time; it never reaches the caller under default
#                  arguments (checked on Chronos 2.3.15).
# 5. Copy-number profile: cell lines x genes, applied AFTER training (see below).

chronos.check_inputs(                      # fail fast on schema and orientation
    readcounts={'screen': readcounts_df},  # rows = sequence_ID, columns = sgRNA
    guide_gene_map={'screen': guide_gene_map},
    sequence_map={'screen': sequence_map},
)

model = chronos.Chronos(
    sequence_map={'screen': sequence_map},
    guide_gene_map={'screen': guide_gene_map},
    readcounts={'screen': readcounts_df},
    negative_control_sgrnas={'screen': negative_control_sgrnas},
)
model.train(nepochs=301)
gene_effects = model.gene_effect                      # attribute, not a method call; lines x genes

# Copy-number correction is a separate post-hoc step, not a constructor argument. It needs at
# least 3 cell lines (RuntimeError below that), the CN frame must be lines x genes and cover
# every gene in gene_effects, and it returns two objects: the corrected matrix and the
# per-gene CN shifts it fitted.
gene_effects_cn, cn_shifts = chronos.alternate_CN(gene_effects, copy_number_df)
gene_probabilities = get_probability_dependent(gene_effects_cn, negative_control_genes, positive_control_genes)
```

**DepMap convention:** A gene-effect score <-1 corresponds to "essential" in that cell line; <-0.5 is "depleting." Each DepMap release (quarterly) provides Chronos gene effects and probabilities.

**Critical:** Chronos benefits most from longitudinal data (multiple timepoints per cell line) but can run with multiple cell lines at a single timepoint, and it trains on a single line with replicates too. Copy number is optional to *training*: Chronos trains without it and `alternate_CN` applies the correction afterwards. The correction itself has a hard minimum of 3 cell lines, so for one or two lines, CN-correct with CRISPRcleanR (pre-hoc, on the counts) and use Chronos only for the gene-effect scores, or not at all.

## CERES (Legacy, Superseded by Chronos)

**Meyers RM et al 2017 *Nat Genet* 49:1779** introduced the first formal CN-correction method at DepMap scale. CERES decomposes per-sgRNA LFC as `sgRNA_efficacy * gene_effect - CN_term(copy_number)`, fitting jointly. Superseded by Chronos at DepMap in 2021 due to Chronos' better handling of screen quality and longitudinal data. CERES remains useful for cross-validation.

## Detect Uncorrected CN Bias

**Goal:** Verify that copy-number bias is corrected (or detect it in raw data).

**Approach:** For genes with matched CN profile, compute Spearman ρ between gene-level LFC and copy number. A negative correlation (-ρ) indicates amplified genes are depleted, i.e., CN artifact.

```python
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

def detect_cn_bias(gene_lfc_df, cn_df, amplified_cn=4, diploid_cn=(1.5, 2.5), low_power_n=8):
    '''Test whether amplified genes are depleted relative to diploid ones.

    The genome-wide Spearman rho alone is not enough: a real focal amplicon covers a handful of
    genes out of ~18,000, so rho stays near zero while those genes are severely depleted. On a
    real HAP1 screen with a planted 8-gene 17q12 amplicon at CN 15, rho was -0.034 (not
    "biased") while amplified genes averaged LFC -2.64 against -0.18 for diploid genes. Report
    the stratified gap, and treat either signal as bias.

    n_amplified>=3 is only the floor for the Mann-Whitney test to run at all, not a guarantee it
    has power to detect a real effect. On the real HT-29 FAM84B/MYC/POU5F1B block (CN=8)
    post-CRISPRcleanR, the residual gap was still -0.98 logFC (MYC itself was left uncorrected;
    its two flanking genes were), but n=3 gave p=0.208 -- bias_present reads False even though
    the bias is real and substantial. Bootstrapping that same real 3-gene mixture at larger n
    shows why: P(p<0.01) is only ~0.04 at n=3, ~0.13 at n=8, ~0.30 at n=15 -- the test needs
    dozens of amplified genes at this effect size before it reliably clears its own threshold.
    Below `low_power_n` amplified genes, do not read bias_present=False as "correction
    succeeded" -- check `suspicious_despite_ns` and the raw `amplified_vs_diploid_gap` instead.
    '''
    merged = gene_lfc_df.merge(cn_df, on='gene')
    rho, p = spearmanr(merged['copy_number'], merged['lfc'])
    amplified = merged[merged['copy_number'] > amplified_cn]['lfc']
    diploid = merged[merged['copy_number'].between(*diploid_cn)]['lfc']
    gap, p_gap = float('nan'), float('nan')
    if len(amplified) >= 3 and len(diploid) >= 3:
        gap = amplified.mean() - diploid.mean()
        p_gap = mannwhitneyu(amplified, diploid, alternative='less').pvalue
    low_power = len(amplified) < low_power_n
    return {
        'cn_lfc_rho': rho,
        'p_value': p,
        'n_amplified_genes': len(amplified),
        'amplified_mean_lfc': amplified.mean(),
        'diploid_mean_lfc': diploid.mean(),
        'amplified_vs_diploid_gap': gap,          # negative = amplified genes more depleted
        'p_amplified_more_depleted': p_gap,
        'bias_present': bool((rho < -0.1 and p < 0.01) or (gap < -0.5 and p_gap < 0.01)),
        'low_power_floor': bool(low_power),        # test ran, but underpowered below low_power_n
        'suspicious_despite_ns': bool(low_power and not pd.isna(gap) and gap < -0.5),
    }
```

Run it once genome-wide and once per candidate amplicon (pass only that region's genes plus the
diploid background), because a single amplicon is invisible in the genome-wide statistic.

**Power caveat:** most candidate amplicons have only 3-10 member genes, and `n>=3` is the floor
for the focal Mann-Whitney test to run at all -- not a guarantee it can detect a real effect at
that size (`low_power_floor: true` in the output below `low_power_n=8`). Re-running this exact
function on real post-CRISPRcleanR HT-29 data at a true 3-gene amplicon block gave
`bias_present: False` (p=0.208) for a residual gap of -0.98 logFC units that was, in fact, still
substantial. Below `low_power_n` amplified genes, do not treat `bias_present: False` as proof
the correction worked; check `suspicious_despite_ns` and the raw `amplified_vs_diploid_gap`
instead, and treat anything below -0.5 as still suspicious even when `p_amplified_more_depleted`
is not significant.

**Threshold (operational convention):** Spearman ρ <-0.10 between LFC and CN indicates genome-wide CN bias, and an amplified-vs-diploid LFC gap below -0.5 with a significant one-sided test indicates focal bias the correlation misses. Either one means correct before hit calling. Run this diagnostic before AND after correction.

## Reconciliation: When CN Correction Fails

If post-CRISPRcleanR or post-Chronos the CN-LFC Spearman is still significantly negative, the correction is incomplete. Possible causes:

1. **Insufficient CN resolution:** A specific 4-copy region went undetected. Refine CN profile with deeper WGS.
2. **CRISPRcleanR position-based correction missed it:** The amplification is small relative to the segmentation algorithm's resolution. Use Chronos with matched CN profile.
3. **Genomic rearrangement creates a "ghost" amplification:** A complex rearrangement appears as normal CN but Cas9 cuts at multiple sites due to translocation breakpoints. Combine WGS structural variants with the analysis.
4. **Cell line has an unusually strong cut-toxicity response:** The artifact may persist; use CRISPRi screens for that line.
5. **`alternate_CN` fits one CN-effect curve across the whole panel:** it is a fitted average, not a per-line correction, so a single extreme-CN line can retain residual signal even after correction. On a real 3-line CN dose-response panel (CN=2/8/15), the CN=15 line's amplicon genes moved from -1.92 to only -1.20 -- still essential-looking -- while the essentials stayed untouched. Re-run `detect_cn_bias` per cell line after `alternate_CN`; compare those results with a pooled summary, but do not assume either view will always conceal the other.

## Apply CN Correction to Pipeline

**Workflow:**

```
1. mageck count (raw counts)
2. screen-qc verification
3. CN diagnostic: Spearman of LFC vs CN (if CN profile available)
4. If bias detected:
   a. CRISPRcleanR (pre-hoc) -> corrected counts -> MAGeCK / BAGEL2 / drugZ
   OR
   b. Chronos (joint model with CN profile) -> gene effects directly
5. Re-diagnose: Spearman of CORRECTED LFC vs CN should be near zero
6. Hit calling
```

For DepMap-style large panels:
```
Chronos handles batch + CN + screen quality in one step; no pre-correction needed.
```

For Project Score-style panel (Behan 2019):
```
CRISPRcleanR was used historically; cross-check with Chronos when CN profile available.
```

## Failure Modes

### CRISPRcleanR removes legitimate essential signal

**Trigger:** A genuine essential gene happens to lie in a region with adjacent uncorrected non-essential signal; the segment-based correction includes the essential.
**Mechanism:** CRISPRcleanR's `ccr.GWclean()` segments sgRNAs by position; segments containing multiple genes with directional consistency are corrected as a unit.
**Symptom:** A known essential drops out of post-correction hit list.
**Fix:** Inspect segments manually; if a known essential was within a corrected segment, investigate. Cross-check with non-CN-corrected MAGeCK + BAGEL2 to see if essential was a hit pre-correction.

### Chronos refuses to CN-correct fewer than 3 cell lines

**Trigger:** One or two cell lines, with or without a matched CN profile.
**Mechanism:** `alternate_CN` fits the CN-to-gene-effect curve across cell lines, so it needs a panel; its own code raises below 3.
**Symptom:** `RuntimeError: Correct for CN should not be used with fewer than 3 cell lines. Consider preprocessing with CRISPRCleanR` -- raised after a full training run, so budget for it before you spend the epochs.
**Fix:** CN-correct the counts with CRISPRcleanR pre-hoc, then hit-call; Chronos training itself works on a single line (with `negative_control_sgrnas`) and its uncorrected gene effects are still useful for QC.

### Spearman ρ still negative after CRISPRcleanR

**Trigger:** Amplification is too small or complex for the segment-based approach.
**Mechanism:** CRISPRcleanR detects systematic spatial patterns; isolated 4-copy regions can slip through.
**Symptom:** Post-correction Spearman ρ -0.05 to -0.10 between LFC and CN.
**Fix:** Refine CN profile (deeper WGS); apply Chronos with matched CN as alternative; or supplement with focal-amplification-aware methods.

### Cell line lacks matched CN profile

**Trigger:** Newly characterized line or rare patient-derived line; WGS not done.
**Mechanism:** CRISPRcleanR needs no CN profile; Chronos trains without one too, but its `alternate_CN` correction cannot run without it.
**Symptom:** Chronos gene effects still carry the artifact; CRISPRcleanR corrects unsupervised.
**Fix:** Run SNP-array (cheap, fast) or low-coverage WGS to obtain CN profile; in interim, use CRISPRcleanR unsupervised mode.

### CN amplification at non-coding region drives apparent essentiality

**Trigger:** Amplification at a gene-poor region; sgRNAs at edge genes get artifactually depleted.
**Mechanism:** Even non-essential genes adjacent to amplifications are depleted because the Cas9 cuts are at the amplified loci.
**Symptom:** Non-essential genes near amplification show LFC <0.
**Fix:** Inspect chromosomal position of "essential" hits; flag genes within 100 kb of known amplifications for orthogonal validation. This is the classic Aguirre 2016 observation.

## CRISPRi/a Alternative

**For variant-function or non-cancer-line essentiality screens, switching to CRISPRi (catalytically dead dCas9-KRAB) avoids the artifact entirely.** No DNA double-strand breaks = no DNA-damage G2 arrest = no copy-number-driven depletion.

| Approach | CN artifact | When to use |
|----------|--------------|-------------|
| Cas9 KO | YES; requires correction | Loss-of-function essentiality, traditional screens |
| CRISPRi | NO | Cancer lines with focal amps; knockdown of cuttable-toxic genes |
| CRISPRa | NO | Gain-of-function; activation screens |
| Base editing | Reduced (single-strand nick) | Variant function |
| Prime editing | Reduced | Precise edits |

See [[library-design]] for CRISPRi (Dolcetto) and CRISPRa (Calabrese) library options.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Spearman ρ (CN vs LFC) | <-0.10 -> bias present | Operational convention |
| Copies for detectable artifact | >6 | Operational convention; response scales with copy number (Aguirre 2016) |
| CRISPRcleanR `min_reads` | 30 (default) | Iorio 2018; lower thresholds in low-coverage screens |
| Chronos gene-effect threshold for "essential" | <-1 (cancer line) | DepMap convention |
| Chronos gene-probability for "essential" | >0.5 | DepMap convention (dependency-probability cutoff) |
| Post-correction Spearman ρ | abs(ρ) <0.05 | Acceptable correction quality |
| Cell-line CN profile resolution | ≥SNP-array level | Below this, CRISPRcleanR unsupervised |
| Cell lines needed for `chronos.alternate_CN` | ≥3 | Chronos 2.3.15 raises RuntimeError below this |
| Amplified-vs-diploid LFC gap | < -0.5 -> focal bias present | Catches single amplicons the genome-wide rho misses |
| `detect_cn_bias` low-power floor | <8 amplified genes -> distrust `bias_present: False` | Real HT-29 3-gene block: p=0.208 for a genuine -0.98 gap |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `RuntimeError: Correct for CN should not be used with fewer than 3 cell lines` | `alternate_CN` needs a panel | CN-correct with CRISPRcleanR instead; Chronos can still score the screen |
| `AssertionError: ... Is your data transposed?` | readcounts passed as guides x samples | Transpose: rows = sequence_ID, columns = sgRNA |
| `ValueError: excess_variance was passed as dict without key for '<library>': {}` at `chronos.Chronos(...)` construction (not `train()`) | `negative_control_sgrnas` not supplied | Pass a dict of non-targeting/non-essential sgRNAs per library |
| `ValueError: set of negative_control_sgrnas is empty` at `chronos.Chronos(...)` construction | `negative_control_sgrnas` was supplied, but one library's list is empty | Populate every library with real non-targeting/non-essential sgRNAs; do not use an empty placeholder list |
| CRISPRcleanR removes a known essential | Segment-based over-correction | Manually inspect segments; cross-check with non-corrected |
| Spearman ρ still -0.15 after correction | Method too coarse for the amp | Refine CN profile; use Chronos |
| ERBB2 listed as essential in SK-BR-3 | Uncorrected HER2 amplification | Always apply correction before hit calling |
| CN profile missing for newly characterized line | Profile not generated | Run SNP-array / low-coverage WGS |
| Hits restricted to non-amplified regions only | Over-correction | Reduce CRISPRcleanR aggressiveness; check known biology |

## References

- Aguirre AJ et al. 2016. *Cancer Discov* 6:914. Copy-number gene-independent toxicity.
- Munoz DM et al. 2016. *Cancer Discov* 6:900. CN amplification CRISPR artifacts.
- Haapaniemi E et al. 2018. *Nat Med* 24:927. Cas9 cutting induces a p53-mediated DNA-damage response.
- Ihry RJ et al. 2018. *Nat Med* 24:939. p53 inhibits Cas9 engineering in human pluripotent stem cells.
- Meyers RM et al. 2017. *Nat Genet* 49:1779. CERES; first formal CN correction at DepMap scale.
- Iorio F et al. 2018. *BMC Genomics* 19:604. CRISPRcleanR.
- Dempster JM et al. 2021. *Genome Biol* 22:343. Chronos.
- Behan FM et al. 2019. *Nature* 568:511. Project Score with CRISPRcleanR-corrected data.
- Pacini C et al. 2021. *Nat Commun* 12:1661. Integrated cross-study dependencies; DepMap quality scoring.
- DepMap Q4 2024+ data releases. https://depmap.org/portal/

## Related Skills

- crispr-screens/screen-qc - CN-LFC Spearman diagnostic; pre-correction QC
- crispr-screens/library-design - Switch to Dolcetto (CRISPRi) to bypass artifact
- crispr-screens/mageck-analysis - MAGeCK on CRISPRcleanR-corrected counts
- crispr-screens/bagel-essentiality - BAGEL2 on CRISPRcleanR-corrected counts
- crispr-screens/hit-calling - Cancer-line hit calling with Chronos
- crispr-screens/batch-correction - Chronos handles batch + CN jointly
- crispr-screens/jacks-analysis - JACKS does not handle CN bias
- clinical-databases/clinvar-lookup - Variant annotation downstream
- copy-number/copy-ratio-segmentation - CN profile derivation upstream
