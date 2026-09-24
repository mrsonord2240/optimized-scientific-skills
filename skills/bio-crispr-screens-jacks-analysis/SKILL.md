---
name: bio-crispr-screens-jacks-analysis
category: Data Analysis
description: Runs JACKS (Joint Analysis of CRISPR/Cas9 Knockout Screens; Allen et al 2019 Genome Research) which models per-sgRNA log-fold-change as the product of a treatment-dependent gene-essentiality term and a treatment-independent guide-efficacy term. Covers the Bayesian decomposition math, the hierarchical efficacy prior shared across screens performed with the same library, when JACKS outperforms MAGeCK (multi-screen joint analysis, libraries with broad efficacy variance) and when it does not (single screen, novel libraries with no prior efficacy), library-reuse efficacy transfer, downstream essentiality interpretation, and the 2.5x sample-size reduction enabled by efficacy-aware testing. Use when running multiple screens with the same library, when guide-level noise is suspected to dominate per-gene signal, when reusing published essentiality reference screens for efficacy priors, or when comparing screens performed across cell lines that share library but differ biologically.
tool_type: python
primary_tool: JACKS
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: JACKS 0.2 (felicityallen/JACKS; checked 2026-09-16), pandas 2.2+, numpy 1.26+, scipy 1.12+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `python run_JACKS.py --help` (run from `JACKS/jacks/` after clone)
- Python: from jacks.jacks_io import runJACKS; help(runJACKS)
- GitHub: install via `git clone https://github.com/felicityallen/JACKS && cd JACKS/jacks && pip install .` (setup.py is in `jacks/`, not the repo root; the PyPI `jacks` is an unrelated package)

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## JACKS CRISPR Screen Analysis

**"Analyze CRISPR screens with guide-level efficacy modeling"** -> Jointly model per-sgRNA log-fold-change across one or more screens as the product of gene essentiality and guide efficacy, sharing efficacy across screens with the same library so that low-quality guides are down-weighted automatically.

- CLI: `python run_JACKS.py countfile replicatefile guidemappingfile [options]` (run from `JACKS/jacks/` after clone)
- Python: `from jacks.jacks_io import runJACKS` for programmatic use; lower-level `from jacks.infer import inferJACKS`
- Output: per-gene effect (`Gene` + one column per cell line) with a matching posterior-std file, per-sgRNA efficacy (`sgrna`, `X1`, `X2`); p-values only with `--ctrl_genes` (see below)

## The JACKS Model (under the hood)

**Why this matters for postdoc-level use:** JACKS decomposes the observed per-sgRNA log-fold-change as:

```
LFC[i, c] = gene_effect[g(i), c] * guide_efficacy[i] + noise
```

where `i` is sgRNA index, `c` is screen condition, `g(i)` is the gene targeted by sgRNA i. Gene effect varies by condition (different cell lines, different treatments) but guide efficacy is intrinsic to the sgRNA sequence and is treated as constant across screens. The model fits both parameters via variational Bayes with hierarchical priors:

- `guide_efficacy[i] ~ Normal(1, 1)` (Gaussian prior, mean 1, variance 1; `mu0_x`, `var0_x` in `jacks.infer.inferJACKSGene`), shared across all sgRNAs
- `gene_effect[g, c] ~ Normal(0, 1e4)` per condition (effectively flat; `mu0_w`, `var0_w`), unless `apply_w_hp` replaces it with a prior fitted to that gene's effects across conditions
- a per-observation noise precision `tau` whose prior weight is `tau_prior_strength=0.5`

The variational posterior gives expected guide efficacy and gene effect. Gene-level p-values are not a likelihood-ratio test: they compare each gene's effect with pseudo-genes built by resampling the guides of supplied negative-control genes (`--ctrl_genes` plus `n_pseudo` > 0).

Initialisation is fixed (efficacy = 1, gene effect = median LFC) and nothing in the inference is random, so gene effects and efficacies are identical across reruns on the same input. Only the pseudo-gene p-values are random (Python's `random`, drawing from the control genes held in a `set` of strings). To make them reproducible, set both `random.seed(<int>)` before `runJACKS` and the environment variable `PYTHONHASHSEED=<int>` before Python starts (checked on JACKS 0.2: `random.seed` alone still gave different p-value files across separate processes; with both fixed two runs matched byte for byte, and a different `PYTHONHASHSEED` changed the file).

**Critical assumption:** Guide efficacy is treated as cell-line independent within the same chemistry. Allen 2019 reports per-sgRNA Cas9 KO efficacy is consistent across randomly selected batches of cell lines (within-chemistry), supporting library-shared efficacy. **However**, efficacy is NOT shareable across chemistries: Cas9 KO efficacy != CRISPRi knockdown efficiency != CRISPRa activation efficiency. JACKS must be run separately per chemistry; use only within the same chemistry on the same library.

## When JACKS Outperforms MAGeCK and BAGEL2

| Scenario | Advantage | Expected gain (Allen 2019) |
|----------|-----------|------------------------------|
| Multi-screen joint analysis (>=3 screens with same library) | Efficacy shared; noise averaged | ~21% lower error vs MAGeCK; 9% vs original BAGEL (Allen 2019 did not benchmark BAGEL2); 91-99% of cell lines improved (method-dependent) |
| Reusing public reference screens (DepMap, Project Score) as efficacy prior (`references/efficacy-prior-and-diagnostics.md`) | Transfer learning | New screens can be smaller; efficacy priors transfer across same-library screens |
| Libraries with broad efficacy variance (e.g. older GeCKOv2) | Down-weights known weak guides | Larger gain than on Brunello (already efficacy-filtered) |
| Heterogeneous quality (mixed plasmid quality across screens) | Per-screen noise estimation | Cleaner per-condition gene effects |

## When JACKS Is Not the Right Tool

- **Single screen, no prior efficacy:** JACKS has nothing to leverage; MAGeCK or BAGEL2 work as well.
- **Single timepoint / two-condition essentiality:** RRA or BAGEL2 simpler and equivalent.
- **Heavy-selection drug screens:** drugZ explicit for chemogenomic; JACKS less sensitive.
- **Cancer-cell-line copy-number screens:** Chronos preferred; jointly models CN bias + screen quality; JACKS does neither.
- **Heavy selection (>40% of guides change):** RRA fails (use MAGeCK MLE); BAGEL2 stays robust.

## Run JACKS Joint Analysis

**Goal:** Jointly analyze multiple CRISPR screens performed with the same library and chemistry.

**Approach:** Provide a count matrix with all samples across all screens, a replicate map identifying which samples belong to which screen and condition, and a sgRNA-to-gene map. JACKS learns guide efficacy shared across screens and gene effects per screen.

Input files (all tab-separated):

- `counts.txt`: rows = sgRNA; first column `sgRNA` (or custom, see `--sgrna_hdr`), then one column per sample.
- `replicatemap.txt`: WITH header `Replicate`, `Sample`, `Control`; column names match the header flags below.
- `guidemap.txt`: WITH header `sgRNA`, `Gene` (the count matrix itself works if it has both columns; point `sgrna_hdr`/`gene_hdr` at them).

```
Replicate                Sample          Control
Screen1_T1               Screen1_T       Screen1_C
Screen1_T2               Screen1_T       Screen1_C
Screen1_C1               Screen1_C       Screen1_C
Screen2_T1               Screen2_T       Screen2_C
Screen2_T2               Screen2_T       Screen2_C
Screen2_C1               Screen2_C       Screen2_C
```

```bash
# Programmatic run through jacks.jacks_io.runJACKS (apply_w_hp stays off; see the variant below)
python scripts/run_jacks_joint.py counts.txt replicatemap.txt guidemap.txt --outprefix jacks_out
# p-value file: add --ctrl-genes NEGv1.txt --n-pseudo 2000 (the Python API defaults n_pseudo=0, the CLI 2000;
# with n_pseudo=0 no p-value file is written). Reproducible p-values: PYTHONHASHSEED=1 ... --seed 1
```

**Variant: hierarchical gene-effect prior (deliberate use only).** `apply_w_hp=True` (CLI `--apply_w_hp`, script `--apply-w-hp`) re-fits the gene-effect prior to each gene's effects across conditions, shrinking them towards each other. The JACKS help marks it "not recommended, use with caution", and it changes rankings materially: on JACKS' own 13-cell-line Project Score example, per-line Spearman rho between the two settings was 0.75-0.89. Use it only when you intend cross-condition shrinkage, and report which setting you used.

```bash
# Equivalent CLI run (run from JACKS/jacks/ after clone; --apply_w_hp stays off)
# --ctrl_sample_hdr names a per-sample control column; use --common_ctrl_sample <name> for one shared control
python run_JACKS.py \
    counts.txt \
    replicatemap.txt \
    guidemap.txt \
    --rep_hdr Replicate \
    --sample_hdr Sample \
    --ctrl_sample_hdr Control \
    --sgrna_hdr sgRNA \
    --gene_hdr Gene \
    --outprefix jacks_out
# Outputs:
#   jacks_out_gene_JACKS_results.txt      gene effect: header `Gene` + one column per cell line
#   jacks_out_gene_std_JACKS_results.txt  matching posterior std per gene per cell line
#   jacks_out_gene_pval_JACKS_results.txt p-values (only with --ctrl_genes; the CLI's --n_pseudo defaults to 2000)
#   jacks_out_grna_JACKS_results.txt      sgRNA-level: header `sgrna`, `X1`, `X2`
#   jacks_out_JACKS_results_full.pickle  full posterior for downstream
#   jacks_out_logfoldchange_means.txt, _logfoldchange_std.txt  raw per-sgRNA log-fold-change vs
#                                          control, mean and std across replicates; written unless --reffile is given
```

## Output Interpretation

| Column | Meaning | Direction |
|--------|---------|-----------|
| `X1` (gene file) | Posterior mean of gene_effect | Negative = essential (depleted); positive = enriched |
| gene std file | Posterior std of the gene effect | Lower = more confident; combine as effect/std for a z-like statistic |
| `X1` (sgRNA file) | Posterior mean of guide efficacy | Centred near 1 and unbounded; the reference Avana set spans negative values to >100 |
| `X2` (sgRNA file) | Second moment E(X^2) of efficacy; std = sqrt(X2 - X1^2) | Confidence in the efficacy estimate |

**Interpretation rule:** A gene is essential if its effect is negative and large relative to its posterior std: divide each cell-line column of the gene file by the same column of the gene std file; effect/std < -2 is roughly a 95% credible deviation from zero. The gene file has no `X1`/`X2` columns. Those belong to the sgRNA file, where `X2` is a second moment, not a std. Supply `--ctrl_genes` (with `n_pseudo` > 0 in Python) to also get a p-value file. Sort by effect (most negative first) for essentiality rank.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Hit call | gene effect < 0 AND abs(effect/std) > 2, from the gene and gene std files | Bayesian z-equivalent (~95% credible); p-values need --ctrl_genes |
| Low-efficacy guide flag | X1 (sgRNA) <0.3 | Operational convention; below this, guide likely non-functional |
| Reference for prior reuse | DepMap or Project Score panel (same library and chemistry only) | Established efficacy distribution |
| Minimum screens for joint efficacy benefit | 3+ | Below this, single-screen tools (MAGeCK/BAGEL2) equivalent |
| Iterations for variational inference | at most 50 per gene, early stop at lower-bound change < 0.1 (JACKS 0.2) | Not exposed on the CLI; see convergence failure mode in `references/failure-modes.md` |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Genes missing from output | sgRNA-to-gene map mismatch (unmatched guides are dropped) | Verify naming consistency; check `len(gene_results) == n_genes_expected` |
| Median efficacy <0.2 | Wrong chemistry assumed by prior | Use a matched-chemistry `--reffile` or override the efficacy prior (see `references/failure-modes.md`) |
| Lower bound still changing at `Iter 50/50` in the DEBUG log | Iteration cap reached | Refit with a higher `n_iter` (see `references/failure-modes.md`) and compare effects |
| p-values differ between runs | Pseudo-genes are sampled with Python `random` from a set of control genes | `random.seed(<int>)` before `runJACKS` AND `PYTHONHASHSEED=<int>` at interpreter start; gene effects themselves are deterministic |
| `<sgRNA> has no sgrna reference in <reffile>` | `--reffile` from a different library | Match library exactly |

## Reference Files

- `references/efficacy-prior-and-diagnostics.md` - read to build a `--reffile` efficacy prior from a reference panel, or to flag low-efficacy guides and genes for library re-design (`scripts/` holds the runnable helpers).
- `references/failure-modes.md` - read when efficacy collapses, results look noisy across cell lines, genes hit the iteration cap, genes go missing, or a `--reffile` raises or seems wrong; holds the CRISPRi hyperparameter and `n_iter` overrides.
- `references/tool-comparison.md` - read to compare JACKS with MAGeCK/BAGEL2 or to reconcile hits where the tools disagree.
- `scripts/run_jacks_joint.py` (Python-API joint run) and `scripts/efficacy_summary.py` (low-efficacy guide/gene summary); `examples/run_jacks.py` wraps the CLI run and the result analysis and plots.

## References

- Allen F et al. 2019. *Genome Research* 29:464. JACKS; original Bayesian joint analysis paper.
- Allen F, Parts L (Wellcome Sanger Institute). https://github.com/felicityallen/JACKS. Official repository.
- Behan FM et al. 2019. *Nature* 568:511. Project Score CRISPR panel; library-wide reference screen data.
- Meyers RM et al. 2017. *Nat Genet* 49:1779. Avana CRISPR DepMap; reference panel for efficacy transfer.

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK RRA/MLE comparison
- crispr-screens/bagel-essentiality - Alternative for essentiality without efficacy modeling
- crispr-screens/library-design - sgRNA design rules informed by JACKS efficacy output
- crispr-screens/copy-number-correction - Chronos preferred for cancer-line multi-screen analyses
- crispr-screens/screen-qc - Pre-JACKS QC; replicate Pearson must pass before joint analysis
- crispr-screens/hit-calling - Cross-method decision tree
- crispr-screens/batch-correction - JACKS does not adjust for batch; pre-correct if necessary
