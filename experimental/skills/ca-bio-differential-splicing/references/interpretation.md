## Per-Tool Failure Modes

### leafcutter: Cluster Mis-Topology

**Trigger:** A cluster spans a complex topology (cassette + alternative donor in same cluster).

**Mechanism:** Cluster-level p-value reports "something in this cluster differs" but doesn't indicate which intron drove the change; downstream analysis needs per-intron effect sizes.

**Symptom:** Significant cluster, multiple introns with different ΔPSI directions, ambiguous biological interpretation.

**Fix:** Inspect cluster in leafviz; report per-intron effect sizes from `ds_results_effect_sizes.txt`; map cluster topology to canonical SE/A5SS/A3SS via flanking exon coordinates manually.

## Reconciliation: When Tools Disagree

The two most common short-read tools answer slightly different questions: rMATS classifies on annotated event templates; leafcutter classifies on observed cluster usage. Disagreement is informative.

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| rMATS sig, leafcutter not sig | rMATS junction imbalance OR rMATS event hits annotation that leafcutter clustered differently | Inspect locus in IGV; check Shiba on the same locus |
| leafcutter sig, rMATS not sig | Novel junction not in rMATS annotation; rMATS `--novelSS` may have missed it | Verify `--novelSS` was on; rerun if not |
| Both sig, opposite ΔPSI direction | First check sign conventions: rMATS `IncLevelDifference` = group 1 - group 2; SUPPA2 dPSI and leafcutter `deltapsi` = group 2 - group 1; Shiba dPSI = alternative - reference group (verified on planted data). Then event class mismatch (rMATS calls SE positive, leafcutter sees A5SS shift in same cluster) | Flip one sign; if still opposite, manually map cluster topology to event class |
| Both sig, same direction | High-confidence call | Report; cross-validate with sashimi-plot |
| All tools null but biology suggests change | Underpowered design or wrong regime | Increase replicates; check whether outlier-splicing-detection regime applies |

**Operational rule:** for high-confidence reporting, require concordant detection in two tools from different algorithmic families (event-based + cluster-based, or LSV + isoform-based). Document both calls and any explainable disagreements.

## rMATS Output Columns Reference

| Column | Meaning |
|--------|---------|
| IJC_SAMPLE_1 / SJC_SAMPLE_1 | Comma-delimited inclusion / skipping junction counts per replicate, group 1 |
| IJC_SAMPLE_2 / SJC_SAMPLE_2 | Same for group 2 |
| IncFormLen / SkipFormLen | Effective lengths normalizing PSI for differential mapping opportunity |
| upstreamES/EE, downstreamES/EE | Flanking exon coordinates (genomic order; strand-agnostic in column meaning) |
| exonStart_0base / exonEnd | Cassette exon coordinates (0-based half-open) |
| PValue | LRT p-value of \|ΔPSI\| > cutoff |
| FDR | BH-adjusted PValue within event class |
| IncLevel1, IncLevel2 | Comma-delimited per-replicate PSI values |
| IncLevelDifference | mean(IncLevel1) - mean(IncLevel2); sign matches --b1 - --b2 order |

## Replicate Count and Power

| Design | Expected power for ΔPSI=0.2 |
|--------|------------------------------|
| n=2 vs n=2 | Marginal; many real effects missed; SUPPA2 unusable |
| n=3 vs n=3 | Adequate at moderate coverage; standard |
| n=5 vs n=5 | Good; recommended for publication |
| n=10+ vs n=10+ heterogeneous | MAJIQ-HET designed for this scale |
| Single patient vs n=20+ controls | Outlier regime (leafcutterMD, FRASER2); see outlier-splicing-detection |

For an effect-size of |ΔPSI|=0.10 (typical biological signal), power generally requires n>=4 and >=20 junction reads per replicate. Below this, expect to miss most real changes.

## Significance and Effect-Size Thresholds

| Stringency | \|ΔPSI\| | FDR | Use case |
|------------|----------|-----|----------|
| Lenient | > 0.05 | < 0.10 | Discovery, exploratory, hypothesis generation |
| Standard | > 0.10 | < 0.05 | Publication; default reporting threshold |
| Stringent | > 0.20 | < 0.01 | Validation cohort, follow-up targets |

For MAJIQ: posterior probability `P(|ΔPSI| > 0.2) >= 0.95` is roughly equivalent to standard stringency (a different scale from FDR; confirm by simulation when reporting). Always document tool, threshold, and rationale.

**Biologically meaningful ΔPSI varies by context:**
- A poison exon shift of |ΔPSI|=0.10 can halve functional protein (huge biology, modest number).
- A stoichiometric isoform shift of |ΔPSI|=0.10 may be physiologically silent.
- Therapeutic ASO target: SMA nusinersen aims for ΔPSI~+0.30 in SMN2 exon 7.

## Confounder Handling

**Trigger:** sequencing batch, RIN, library prep date or sex correlates with the comparison. rMATS' default LRT takes no covariates (`--paired-stats` handles pairing only), so hits can be driven by batch: PCA on the PSI matrix shows samples clustering by batch rather than group. **Check this first**; test every leading PC against batch and group, not only PC1 (a batch effect can appear on PC2 or later).

**Batch identical to (or perfectly nested in) condition cannot be adjusted.** leafcutter_ds.R with `batch = group` exits 0 with no warning but returns non-informative cluster p-values near 1; test the design-matrix rank before trusting any covariate model (the snippet below does).

Workarounds for rMATS:
1. **Stratification**: run rMATS within each batch separately and meta-analyze.
2. **Per-event regression with a group term (logit-transformed PSI)**: PSI is bounded [0,1], so logit-transform, fit `logit_psi ~ group + batch (+ RIN)` and test the **group coefficient**. Do not regress out batch alone and rank-test the residuals: with 3 vs 3 the rank-sum test cannot go below p=0.1 (0/20 strong events reached p<0.05), and residualizing on an imbalanced batch removes group signal.
3. **For n<=4 per group, switch to leafcutter with a batch column first** (R function accepts `confounders` matrix; CLI accepts confounders as additional columns in the groups file). On the checked 4v4 simulation it retained 18/30 true events while removing batch-driven calls; the per-event regression was much less powered.

```python
import numpy as np
import pandas as pd
import patsy
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

# meta: one row per sample (columns group 0/1, batch, RIN); psi: one row per sample for ONE event (columns psi + the same covariates)
design = patsy.dmatrix('group + C(batch) + RIN', meta)
if np.linalg.matrix_rank(np.asarray(design)) < design.shape[1]:
    raise ValueError('a covariate is aliased with group: this comparison cannot be adjusted')

eps = 1e-3  # logit-transform PSI before modelling (PSI is bounded [0,1])
p = psi['psi'].clip(eps, 1 - eps)
psi['logit_psi'] = np.log(p / (1 - p))
fit = smf.ols('logit_psi ~ group + C(batch) + RIN', data=psi).fit()
p_group = fit.pvalues['group']   # one p per event; then multipletests(p_values, method='fdr_bh') across events
```

Power is limited by residual degrees of freedom (samples minus model terms): on a simulation with an imbalanced batch plus RIN, 3 vs 3 found 0/21 strong events at BH q<0.05 (9/21 at nominal p<0.05; residual df=2), 5 vs 5 found 21/21 with 1 false call among 46 null events. With few replicates drop RIN, or use stratification.

**leafcutter** accepts confounders two ways (verified in the script source and by a run):
- **R function**: `differential_splicing(counts, x, confounders=numeric_matrix)` accepts a numeric covariate matrix
- **CLI script**: `leafcutter_ds.R` reads confounders from **additional columns in the groups file** (3rd, 4th, ... columns), NOT from a `--confounders` flag; use labels like `b1`/`b2` for categorical variables (numeric columns are treated as continuous)

**MAJIQ** does not accept arbitrary confounders (its documentation lists a separate `moccasin` batch-correction step; unverified); use stratification or switch tool.

## Multi-Group / Multi-Factor Designs

| Design | Approach |
|--------|----------|
| 3 groups (e.g. drug A, drug B, control) | Pairwise rMATS or leafcutter; OR limma/DESeq2 on logit-PSI matrix |
| Time-course (e.g. 0h, 6h, 24h) | DEXSeq on event counts with time as factor; or limma::lmFit on PSI matrix |
| 2x2 factorial (genotype × treatment) | DEXSeq with interaction term; rMATS pairwise on interaction subsets |
| Continuous covariate (dose, age) | limma::lmFit on logit-PSI ~ covariate |

For complex designs, custom regression on the PSI matrix is more flexible than rMATS/leafcutter pairwise.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `rMATS: numpy.AxisError` | rMATS version mismatch with numpy >=2.0 | Pin numpy<2.0 or update rMATS-turbo to >=4.3 |
| rMATS exits 0, `SE.MATS.JC.txt` has only a header (or no FDR column) | `-t` does not match the reads, `--readLength` mismatch without `--variable-read-length`, or `--paired-stats` without PAIRADISE (`no package called PAIRADISE` in the log) | Fix the flag / install PAIRADISE; assert rows + FDR after every run |
| `leafcutter_ds.R: The number of samples in the smallest group is less than min_samples_per_intron` | Default `-i 5` with a group of <5 samples | `-i N -g N` with N = smaller group size (leafcutter section) |
| `leafcutter_ds.R: undefined columns selected` | Sample names in groups.txt do not equal the counts-table column names (.junc basenames) | Rename groups.txt entries to the .junc basenames |
| leafcutter clustering writes 0 introns | Contig not `chr1..22,X,Y`/`1..22,X,Y`; `-m` too high for the depth; or `.junc` strand `?` (no XS tag) | `-k True`; lower `-m`; use XS-tagged BAMs |
| SUPPA2 `ImportError` on `multipletests` | statsmodels 0.15 | Install statsmodels <=0.14 |
| `SUPPA2: no events with sufficient coverage` | Salmon/kallisto TPM filter too strict upstream | Lower upstream TPM threshold; verify event annotations |
| `regtools: too many open files` | Many BAMs in one batch | `ulimit -n 4096` or batch in groups |

## Result Prioritization

**Goal:** Rank events by combined statistical and biological significance for follow-up.

**Approach:** Composite score combining FDR and effect size, then enrich for biology (RBP binding, NMD sensitivity, conservation, disease relevance).

```python
import pandas as pd
import numpy as np

# `significant` = the filtered rMATS table from the rMATS section
sig = significant.copy()
sig['score'] = -np.log10(sig['FDR']) * sig['IncLevelDifference'].abs()
sig['exon_length'] = sig['exonEnd'] - sig['exonStart_0base']
sig['nmd_likely'] = (sig['exon_length'] % 3 != 0)
top_events = sig.nlargest(50, 'score')
```

Cross-reference top hits with:
- **eCLIP/ENCODE RBP target databases** (POSTAR3, oRNAment, RBP2GO) -> candidate trans-regulators
- **Disease-specific signatures**: SF3B1 cancer signal = cryptic 3'ss ~10-30nt upstream of the canonical one (MDS/CLL/UM); TDP-43 cryptic exons (UNC13A, STMN2, ATG4B) for ALS/FTD
- **Conservation**: VastDB cross-species PSI for evolutionary support
- **Splice-site predictions**: SpliceAI scores for the involved sites (see splice-variant-prediction)

## Common Pitfalls

- **Junction read imbalance** (cassette exon flanks have unequal mapping opportunity) inflates rMATS false positives; Shiba explicitly corrects this.
- **Forgetting NMD direction** — increased PSI of a poison exon decreases protein. Always check whether the alternative form is PTC-introducing using ORF-aware annotation.
- **Cryptic splicing in TDP-43 loss / SF3B1-mutant samples** — annotation-bound tools (rMATS, SUPPA2) miss these; need leafcutter or MAJIQ with denovo mode.
- **Forgetting strand** — wrong `--libType` halves usable junctions. Confirm with RSeQC `infer_experiment.py`.
