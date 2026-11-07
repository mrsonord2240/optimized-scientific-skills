# i2_pseudobulk_de: design notes (not agent-visible)

## Pitfall and source
Multi-donor scRNA-seq DE: treating cells as independent replicates inflates significance; donor is the unit of
replication, so aggregate to donor x condition (pseudobulk) and model with donor as a blocking factor.
Squair et al. 2021 Nat Commun (single-cell DE methods, false discoveries); Murphy & Skene 2022 Nat Commun.

## World (generate.py, seed 1, deterministic)
- 12 donors, 3,367 cells, 4,000 genes, one cell type. Each donor has cells in both conditions; counts per
  donor x condition range 1..340. Composition is deliberately unbalanced: controls come mostly from donors 1-3,7,8, treated mostly from
  donors 4-6 (so donor effects are confounded with condition at the cell level).
- Excluded by the task rule (<5 cells in either condition): 4 donors with (control,treated) = (150,2), (220,3), (4,300), (180,1).
  8 donors are retained; the smallest retained cell is 12 treated cells.
- Model: counts ~ gamma-Poisson, size ~6 per gene, mu = library(lognormal) x fraction_g x exp(donor_dg + sample_dcg + c*lfc_g ln2).
  donor_dg: SD 1.5 ln on a random 20% of genes, 0.6 ln elsewhere. sample_dcg (donor x condition): SD 0.2 ln, all genes.
- Truth = 300 genes with true |log2FC| in [1.4, 2.6], mean expression >= 0.5 UMI/cell. Signs chosen greedily so up/down
  expression mass balances (log2 composition drift -0.0001), otherwise library-size normalisation shifts every null
  gene (first draft: median null LFC -0.74) and the truth becomes normaliser-dependent.
- 150 "modest" genes with true |log2FC| in [0.3, 0.7] (real but below the threshold; not in truth). They catch
  people who drop the fold-change rule (QP paired with p only: 586 genes, precision 0.51).
- truth.json holds the 300 genes, true log2FC, modest genes, and large-donor-effect genes.

## Grader
`marker_gene_precision_recall`, flat list, no expected_count. precision = TP/len(submitted), recall = TP/300.
Thresholds precision >= 0.65, recall >= 0.70 (constructed truth, not solver output).

## Results (seed 1; precision / recall)
| solver | n | prec | rec | outcome |
| --- | --- | --- | --- | --- |
| correct: OLS log2CPM + donor | 263 | 0.985 | 0.863 | pass |
| correct: quasi-Poisson GLM + donor | 333 | 0.895 | 0.993 | pass |
| correct: paired t on log2CPM | 263 | 0.985 | 0.863 | pass |
| cell Wilcoxon (retained donors / all donors) | 700 / 688 | 0.39 / 0.40 | 0.92 / 0.91 | fail (precision) |
| cell Welch t | 712 | 0.386 | 0.917 | fail (precision) |
| pseudobulk unpaired OLS / QP | 137 / 180 | 1.00 / 0.97 | 0.46 / 0.58 | fail (recall) |
| pseudobulk no donor filter, log2CPM OLS or paired t | 174 | 1.00 | 0.58 | fail (recall) |
| pseudobulk no donor filter, quasi-Poisson GLM | 344 | 0.866 | 0.993 | BLIND SPOT: passes |
Seeds 2 and 3 (scratch copies) give the same pattern with similar margins.

Tolerance rationale: correct-variant sweep (prior count 0.25/1/5/20/100 on log2CPM, own vs pooled-count vs mean-of-donor-CPM
LFC, QP) gives precision 0.76-1.0 and recall 0.77-0.99. Traps sit at precision <= 0.40 (cell level) or recall <= 0.58 (unpaired, no-filter).
Thresholds 0.65 / 0.70 fall in the gaps (margins about 0.25 precision, 0.07-0.12 recall). The recall margin on the correct side is the
thinnest: extreme pseudocounts (log1p CP10k-scale, prior 100 CPM) or 0.25 prior reach 0.77-0.80.

## Mechanisms
- Cell-level tests: donor baseline effects (large on 20% of genes) are confounded with condition through the unbalanced
  cell composition, so hundreds of null genes get large pooled fold changes AND tiny p-values. The fold-change rule does not rescue
  them (that is why donor SD on the V genes is large).
- Unpaired pseudobulk: valid but loses power to donor variance, since the donor effect is left in the error (recall about 0.5).
- No-filter pseudobulk on log scale: 2-4-cell pseudobulks (mostly zeros) add huge variance to the paired differences. Count GLMs
  weight those samples by their Poisson variance, so they barely matter (blind spot, honest to record: the filter is
  a variance-control step that only bites for equal-weight log-scale methods).

## Blind spots / ambiguities / gaming
- Blind spot: no-filter pseudobulk with a count GLM (edgeR/DESeq2-like; only the quasi-Poisson proxy was tested; R DESeq2/edgeR/limma are NOT installed
  in either local docker image, so no R variant was run) passes.
- The task defines log2FC as log2 ratio of mean normalised expression (not mean of logs) to remove friction; a mean-of-log1p route
  (scanpy-style per-donor mean of log1p CP10k) gave recall 0.58 in the earlier draft due to compression and would be graded as a fail.
  With a pooled-count LFC and QP, precision drops to 0.76 (donor-composition-biased pooled ratios push modest genes above 1): still passes.
- Paired-vs-unpaired: task states "each donor was sampled in both conditions" (a data fact) but not how to model it.
  A method choosing donor as a random effect (lme4/limma duplicateCorrelation) should work like the paired OLS; untested here.
- Gaming: submitting a huge list fails on precision; submitting only the top ~200 by significance from a pseudobulk fit works
  but is a legitimate high-precision answer. Gene symbols are opaque (SYNGxxxx), no recall from memory possible.
- The 4 excluded donors include one with 300 treated cells and 4 control cells: dropping only "small" donors by total count would be wrong,
  the rule is stated per condition.
- Cell-level with exact scanpy defaults (Wilcoxon, tie_correct false) was approximated by tie-corrected z; no material difference expected.

## Runtime
generate.py ~5 s; solvers.py ~1 min. Requires numpy, scipy, pandas, statsmodels (`pip install statsmodels`).
