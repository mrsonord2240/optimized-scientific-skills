---
name: bio-crispr-screens-screen-qc
description: Quality control for pooled CRISPR screens covering library representation, Gini index, log-skew, replicate Pearson and Spearman concordance, essentialome precision-recall AUC against CEGv2 (Hart 2017), Cas9 cut-toxicity diagnostics, copy-number amplicon detection (Aguirre 2016 / Munoz 2016), bottleneck propagation through plasmid pool, infection, selection, and endpoint stages, MOI verification, and DepMap-style screen-quality scoring. Use when assessing screen quality before hit calling, deciding whether to repeat or rescue a screen, diagnosing low-confidence hits, choosing between MAGeCK / BAGEL2 / Chronos based on quality grade, picking a normalization strategy from QC signatures, or evaluating whether an in-vivo screen retained adequate library complexity.
tool_type: python
primary_tool: MAGeCK-VISPR
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: MAGeCK 0.5+ (count + VISPR), MAGeCKFlute 2.0+ (R), pandas 2.2+, numpy 1.26+, scikit-learn 1.4+, matplotlib 3.8+, seaborn 0.13+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mageck --version` then `mageck count --help`
- R: `packageVersion('MAGeCKFlute')` then `?BatchRemove` / `?FluteRRA`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Install and Inputs

```bash
conda install -c bioconda mageck   # not on PyPI
pip install pandas numpy scipy matplotlib seaborn scikit-learn
# MAGeCK QC dashboard
conda install -c bioconda -c conda-forge mageck-vispr
# R dashboard (optional)
R -e "remotes::install_github('WubingZhang/MAGeCKFlute')"   # removed from Bioconductor at 3.22
```

Required inputs: MAGeCK count output (`screen.count.txt`), plasmid-pool counts (separate file or first sample), known copy-number profile per cell line (from WGS / SNP-array / ASCAT / matched cell-line database), and CEGv2 / NEGv1 reference gene sets (CEGv2 from Hart 2017, NEGv1 from Hart 2014; `hart-lab/bagel` repository).

## CRISPR Screen Quality Control

**"Audit my CRISPR screen quality before hit calling"** -> Assess library representation, replicate concordance, depth, drift, and biological signal recovery using DepMap-grade metrics, then decide whether the screen is usable, salvageable, or must be repeated.

- Python: `pandas` + `scikit-learn` for Gini, AUC, PCA; `MAGeCKFlute` (R) for one-shot QC dashboard
- CLI: `mageck count` writes Gini and mapping stats unconditionally to `<prefix>.countsummary.txt` (`GiniIndex`, `Reads`, `Mapped`, `Percentage`); MAGeCK-VISPR for an interactive dashboard

## QC Stage Hierarchy

A pooled screen has six distinct bottlenecks where complexity can collapse. Audit each:

| Stage | Metric | Acceptable threshold | Failure consequence | Detail |
|-------|--------|----------------------|----------------------|--------|
| Plasmid pool | Gini, skew, % zero-count guides | Gini <0.1, skew <2 (Joung 2017 states <10), zero <0.5% | Missing guides cannot be screened; dropout indistinguishable from non-coverage | Below; [failure modes](references/failure-modes.md) |
| Day-0 infection | Library coverage, MOI verification | ≥99% guide detection at 500x cells/sgRNA; MOI 0.3 | Founder effects; polyclonality with high MOI | [depth-and-moi](references/depth-and-moi.md) |
| Selection (puro/blast) | % cells surviving, time-course Gini | 30-40% survival at 5-7 days; Gini drift <0.05 | Selection artifact; fast-growers enriched | Below |
| Endpoint | Replicate correlation, depth | Pearson >=0.8 on log-counts (MAGeCK-VISPR floor), Spearman >0.7-0.8, >500 reads/sgRNA (Joung 2017 screening) | Noise dominates; FDR inflates | Below; [depth-and-moi](references/depth-and-moi.md); [pca-and-composite-score](references/pca-and-composite-score.md) |
| Biological signal | CEGv2 PR-AUC, NEGv1 false-positive rate | PR-AUC >0.7 at FDR 5% (community "passing" convention); CEGv2 enrichment in top 1k | Screen lacks essentiality signal; hits not credible | [essentialome-recovery](references/essentialome-recovery.md) |
| Copy-number artifact | Amplified-region enrichment, sgRNA-cut-count correlation | No correlation between sgRNA off-target count and depletion | False-positive essentiality at amplicons; ERBB2 in HER2+ etc. | [copy-number-bias](references/copy-number-bias.md) |

Each metric below, and in the reference files, quantifies one of these stages.

## Library Representation Metrics

**Goal:** Detect dropout, oversaturation, and library bottlenecks at each sequencing stage.

**Approach:** Compute per-sample zero-count fraction, low-count fraction (<30 reads, the CRISPRcleanR `ccr.NormfoldChanges` default), and percentile-based skew, then track how these change between plasmid -> Day-0 -> endpoint to localize the bottleneck.

```bash
python scripts/library_representation.py screen.count.txt --out library_representation.tsv
```

`scripts/library_representation.py` prints one row per sample (`n_sgrnas_detected`, `pct_zero`, `pct_lowcount`, `median_count`, `p10_count`, `p90_count`, `skew_ratio`) and also exports `library_representation(counts_df)` and `stage_specific_thresholds()` for import. The stage limits it returns (pct_zero_max / skew_max / gini_max: plasmid 0.5 / 2.0 / 0.10, day_0 1.0 / 2.5 / 0.12, endpoint 5.0 / 10.0 / 0.30) follow Joung 2017 (zero-count, skew) and MAGeCK-VISPR (Gini); skew 2.0 is a stricter modern convention, Joung 2017 states <10.

**Interpretation:** Plasmid pool failing Gini <0.1 indicates synthesis or amplification bias; the screen is unfit for use. Endpoint Gini drifting above 0.30 indicates either heavy biological selection (acceptable for strong-phenotype drug screens) or a bottleneck (must be diagnosed). The Day-0 vs plasmid delta isolates whether the issue arose during infection (cloning is unlikely to lose specific guides between extraction and infection -- the change happens in cells).

## Gini Coefficient

**Goal:** Quantify how unevenly reads are distributed across sgRNAs in a single sample.

**Approach:** Sort non-zero counts ascending, compute Gini via the cumulative-fraction formula. Compare against stage-specific thresholds.

```python
def gini(x):
    '''Gini coefficient: 0 = perfect equality, 1 = maximal inequality.
    Uses non-zero counts only; zero-count sgRNAs handled separately by % zero.'''
    x = np.sort(x[x > 0].astype(float))
    if x.size == 0:
        return np.nan
    n = x.size
    cumx = np.cumsum(x)
    return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
```

**Stage-specific thresholds:** only the plasmid Gini <=0.1 is a published cutoff (MAGeCK-VISPR); the remaining grades are operational convention.

| Stage | Excellent | Acceptable | Concerning | Failure |
|-------|-----------|------------|------------|---------|
| Plasmid pool | <0.10 | <0.15 | 0.15-0.20 | >0.20 |
| Day 0 (post-infection) | <0.12 | <0.18 | 0.18-0.25 | >0.25 |
| Endpoint (post-selection) | <0.30 | <0.40 | 0.40-0.55 | >0.55 |

A Gini that climbs from 0.10 (plasmid) to 0.45 (endpoint) is expected when the screen exerts strong selection (drug, lethal-condition). A Gini that climbs to 0.45 without any biological selection (e.g., a control-vs-control timepoint comparison) indicates technical drift.

## Replicate Concordance

**Goal:** Verify that biological/technical replicates agree before testing for between-condition differences.

**Approach:** Compute pairwise Pearson on log10(counts+1) (MAGeCK-VISPR convention) and Spearman ρ on raw rank, between every replicate pair within a condition. Flag any pair below the MAGeCK-VISPR floor of 0.8 Pearson on log-scale.

```python
def replicate_concordance(counts_df, condition_map):
    '''condition_map: {condition_name: [sample_col1, sample_col2, ...]}.'''
    log_counts = np.log10(counts_df + 1)
    rows = []
    for cond, samples in condition_map.items():
        if len(samples) < 2:
            continue
        for i in range(len(samples)):
            for j in range(i+1, len(samples)):
                r_pearson = log_counts[[samples[i], samples[j]]].corr().iloc[0, 1]
                r_spearman = counts_df[[samples[i], samples[j]]].corr(method='spearman').iloc[0, 1]
                rows.append({'condition': cond, 'rep1': samples[i], 'rep2': samples[j],
                             'pearson_log': r_pearson, 'spearman': r_spearman})
    return pd.DataFrame(rows)
```

**Thresholds:** 0.8 Pearson is the MAGeCK-VISPR floor; the stricter grades are operational convention.

| Metric | Excellent | Acceptable | Failure |
|--------|-----------|------------|---------|
| Pearson on log10(counts+1) | >0.95 | >0.85 | <0.80 |
| Spearman on raw ranks | >0.85 | >0.70 | <0.60 |

**When Pearson is high but Spearman is low**, a few outlier sgRNAs are driving correlation (one extreme guide dominates). Inspect the scatterplot; typically caused by PCR jackpotting at a single guide. Hit calling should use a method that ranks (RRA, drugZ) rather than one that fits per-sgRNA fold change directly.

## Reference Files

Read the one that matches the failing stage or the request; they hold the code and method detail.

| File | Read when |
|------|-----------|
| [references/essentialome-recovery.md](references/essentialome-recovery.md) | Computing CEGv2/NEGv1 PR-AUC, or PR-AUC is low despite good Gini and Pearson |
| [references/copy-number-bias.md](references/copy-number-bias.md) | Cancer-cell-line screen, amplicon genes among the hits, or gating on CN bias (two-rule diagnostic) |
| [references/depth-and-moi.md](references/depth-and-moi.md) | Auditing reads per sgRNA, total-read CV, or verifying MOI and the Poisson multi-guide fraction |
| [references/pca-and-composite-score.md](references/pca-and-composite-score.md) | Checking condition vs batch clustering, or building the single pipeline-gate score |
| [references/failure-modes.md](references/failure-modes.md) | A metric fails and the cause is not obvious (PCR bias, Cas9 heterogeneity, amplicons, outlier replicate, high MOI, CRISPRi TSS) |

Runnable code: `scripts/library_representation.py`, `scripts/essentialome_recovery.py`, `scripts/cn_bias.py` (each takes files as arguments and is importable), and the end-to-end `examples/screen_qc.py` (edit its `stage_map` / `condition_map`, run in the folder holding `screen.count.txt`).

## Interpretation Notes

- Plasmid-pool sequencing (Gini <0.1, >=99% guide detection at >25 reads/guide) is non-negotiable; everything downstream is normalized against this baseline. A screen with an un-sequenced plasmid pool is uninterpretable.
- Day-0 vs plasmid: Pearson >0.9 is expected; below it, diagnose the infection step.
- CEGv2 PR-AUC is the single most diagnostic metric: a pass means the screen has biology even if individual sample metrics look weak, and a screen below 0.5 cannot be fixed in software. "Passing" every earlier stage is necessary but not sufficient; PR-AUC is the final gate.
- Drug screens: endpoint Gini drifting to 0.3-0.5 is normal because biology drives selection. Compare against vehicle, not Day 0, for any chemogenomic interpretation.
- CRISPRi/a: expect lower per-gene PR-AUC than Cas9 (not every essential responds to knockdown as it does to knockout); calibrate against the DepMap CRISPRi sub-essentialome rather than CEGv2.
- Pick the hit-calling method from the quality grade: high quality -> MAGeCK MLE or Chronos; low quality -> RRA or drugZ; cancer line -> Chronos with CN correction; in vivo -> bottleneck-adjusted thresholds.

## When NOT to Use This Skill

Screen QC and the CN-bias diagnostics here assess **research cell-line data**. They do not support a
decision about an individual patient: an amplification that survives correction in a cell line is not
a confirmed biomarker, and no QC metric here speaks to a person's treatment. If a request mixes screen
QC with a clinical decision, answer the QC part and direct the clinical part to validated diagnostics
and a treating physician.

## Input Validation

Before any metric, check the count table and fail with a clear message rather than a traceback:

| Check | Expectation | If it fails |
|-------|-------------|-------------|
| Required columns | an sgRNA identifier column (index) and a `Gene` column, then one numeric column per sample; a default `0..n-1` index means the identifier column was lost | name the missing column; do not guess |
| Dtypes | every sample column numeric | report which column is non-numeric and the first offending value |
| All-zero sample | at least one non-zero count per sample | report the sample as failed at sequencing, and skip (not `nan`-propagate) its Gini and correlation |
| Negative counts | none | reject the file; these are not counts |
| Duplicate sgRNA IDs | none | report the duplicates; MAGeCK-format tables should be unique |

`gini()` above already returns `np.nan` for an all-zero sample rather than raising; keep that guard in
any copy of it.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Plasmid Gini | <0.10 | Li W et al 2015 MAGeCK-VISPR *Genome Biol* 16:281 |
| Plasmid skew ratio (p90/p10) | <10 (Joung 2017); <2 is a stricter modern convention | Joung 2017 *Nat Protoc* 12:828 |
| % zero-count sgRNAs (plasmid) | <0.5% | Joung 2017 *Nat Protoc* 12:828 |
| % zero-count sgRNAs (endpoint) | <1% ideal, <5% tolerated | Li W et al 2015 *Genome Biol* 16:281 |
| Replicate Pearson on log10(counts+1) | >=0.8 (MAGeCK-VISPR floor); >0.95 ideal | Li W et al 2015 MAGeCK-VISPR *Genome Biol* 16:281 |
| Replicate Spearman | >0.70 | Operational convention |
| CEGv2 PR-AUC at FDR 5% | >0.70 passing; >0.85 high quality | Community convention (CEGv2 from Hart 2017 *G3* 7:2719) |
| Reads per sgRNA per sample | >100 plasmid QC and >500 screening (Joung 2017); 300+ (MAGeCK-VISPR) | Joung 2017; Li W et al 2015 |
| Library coverage at infection | 500x cells/sgRNA | Joung 2017; DepMap |
| In-vivo coverage | 50-200x at endpoint | Bottleneck-limited; see [[in-vivo-screens]] |
| MOI at infection | 0.3 strict | Poisson: P(≥2)=4% at 0.3 vs 9% at 0.5 |
| CN-bias Spearman ρ (LFC vs copy number), pre-correction fail | abs(ρ) <0.10 | Operational convention |
| CN-bias Spearman ρ, post-correction target | abs(ρ) <0.05 | Operational convention; the tighter bar after CRISPRcleanR/Chronos |
| CN-bias amplified-vs-diploid LFC gap | > -0.5 | Catches focal amplicons that ρ misses |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Plasmid Gini >0.2 | PCR over-amplification | Re-sequence; cap PCR at 15 cycles |
| Endpoint PR-AUC <0.5 | Cas9 not selected pre-screen | Select Cas9+; redo if mid-screen |
| Pearson high, Spearman low | A few outlier sgRNAs dominate | Use RRA / rank-based hit calling |
| Replicate Pearson <0.8 | Library-prep failure on one rep | Drop outlier; re-run if singleton |
| % zero increases dramatically Day-0 -> endpoint | Selection bottleneck | Reduce selection pressure; or increase coverage |
| Top hits include amplified-region genes | CN bias | CRISPRcleanR or Chronos |
| MOI verification shows 0.6+ | Over-infected | Re-run at lower MOI; no rescue |
| Detected fraction <90% in Day 0 | Coverage too low | Increase cells; expect drift |

## References

- Joung J et al. 2017. *Nat Protoc* 12:828. Genome-wide screen protocol; coverage and depth conventions.
- Li W et al. 2014. *Genome Biol* 15:554. MAGeCK.
- Li W et al. 2015. *Genome Biol* 16:281. MAGeCK-VISPR; Gini, zero-count, depth and replicate-correlation QC cutoffs.
- Wang B et al. 2019. *Nat Protoc* 14:756. MAGeCKFlute; QC dashboard.
- Hart T et al. 2017. *G3* 7:2719. CEGv2 core-essential reference set; PR-AUC screen-quality benchmarking.
- Hart T et al. 2014. *Mol Syst Biol* 10:733. Gold-standard essential and non-essential reference sets; source of NEGv1.
- Aguirre AJ et al. 2016. *Cancer Discov* 6:914. Copy-number amplicon false-essentiality.
- Munoz DM et al. 2016. *Cancer Discov* 6:900. Copy-number gene-independent toxicity.
- Pacini C et al. 2021. *Nat Commun* 12:1661. Integrated cross-study dependencies; NNMD screen-quality metric and cross-study batch correction.
- Meyers RM et al. 2017. *Nat Genet* 49:1779. CERES; mechanism of CN bias.
- Sanson KR et al. 2018. *Nat Commun* 9:5416. Dolcetto/Calabrese TSS rules.
- Dempster JM et al. 2021. *Genome Biol* 22:343. Chronos screen-quality model.

## Related Skills

- crispr-screens/library-design - Compose libraries that pass plasmid QC
- crispr-screens/mageck-analysis - Run MAGeCK count to generate QC inputs
- crispr-screens/copy-number-correction - Remediate Aguirre / Munoz CN artifact
- crispr-screens/batch-correction - Address inter-batch / cell-line confounding
- crispr-screens/hit-calling - Pick method by QC grade
- crispr-screens/in-vivo-screens - In-vivo-specific bottleneck QC
