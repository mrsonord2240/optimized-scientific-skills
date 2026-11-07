# i8_microbiome_da (not agent-visible)

## Pitfall and source
Sequencing returns compositions, not absolute abundances. When a few taxa bloom, the read share of every
other taxon falls, so tests on relative abundance (Wilcoxon / t-test on proportions, rarefied counts, TSS/CPM
+ linear model, library-size-offset NB GLM) call most of the community "differential".
Compositionally aware methods (ANCOM-BC, ALDEx2, LinDA, CLR, reference-frame centring) recover the real changes.
Sources: Gloor et al. 2017 Front Microbiol; Nearing et al. 2022 Nat Commun. Designed from the brief and general
domain practice, no Skill text read before the design was frozen.

## World (generate.py, seed 20260919)
- 200 taxa x 40 samples (20 control / 20 case, order shuffled). Depth log-uniform 8k-40k reads, independent of group.
- ABSOLUTE abundance = lognormal base share (sd 1.6) x group FC x sex effect x per-sample biological lognormal noise
  (sd 0.4-0.8 per taxon, capped 0.6 on truth taxa). Reads = multinomial(depth, composition x Gamma(25) noise).
- Truth (8 taxa): 5 bloomers x10, x8, x7, x6, x5 (combined baseline share 10-20%, drawn from abundance ranks 3-35)
  and 3 decliners x0.25, x0.30, x0.35 (abundant taxa). Every other taxon has identical absolute abundance in both groups.
- Consequence: total load is 2.03x higher in cases; the median unchanged taxon falls to 0.54x relative abundance.
- Sparsity: 6 all-zero taxa, 9 near-zero taxa (a few reads in a handful of samples), ~30% of other null taxa have
  structural zeros (present in 40-90% of samples). 22% of the table is zeros.
- Nuisance: sex (14F/6M in cases, 7F/13M in controls); 10 non-truth taxa carry a sex effect (x2.2-x4). Sex is NOT
  part of the truth. Adjusting is helpful but not required in this realisation.
- Taxon ids are randomly permuted (ASVnnn does not encode abundance or truth). Taxonomy is a random draw of real gut
  lineages; carries no signal. `truth.json` has the answer plus true FCs; never expose it.

## Grader
`marker_gene_precision_recall`, flat list, `answer_field=differential_taxa`, truth = the 8 taxa,
precision_at_k >= 0.40, recall_at_k >= 0.75 (>= 6 of 8). No `expected_count` (variable-length answer).

### Calibration (frozen data)
| method | n called | TP | FP | prec | rec | pass |
| --- | --- | --- | --- | --- | --- | --- |
| CLR (pc 0.5) + OLS(sex), prevalence >= 8 | 9 | 8 | 1 | 0.89 | 1.00 | yes |
| CLR unadjusted | 8 | 8 | 0 | 1.00 | 1.00 | yes |
| CLR + Wilcoxon | 9 | 8 | 1 | 0.89 | 1.00 | yes |
| median-of-log-ratio centring | 10 | 8 | 2 | 0.80 | 1.00 | yes |
| ALDEx2-style (128 Dirichlet MC, mean p) | 7 | 7 | 0 | 1.00 | 0.88 | yes |
| LinDA-like (mode-corrected OLS on log) | 7 | 7 | 0 | 1.00 | 0.88 | yes |
| DESeq-style median-of-ratios size factors + OLS | 10 | 8 | 2 | 0.80 | 1.00 | yes |
| CLR pc=1, NO prevalence filter | 16 | 8 | 8 | 0.50 | 1.00 | yes (near miss) |
| Wilcoxon on proportions (filtered / unfiltered) | 89 | 8 | 81 | 0.09 | 1.00 | fail |
| Welch t on proportions | 89 | 8 | 81 | 0.09 | 1.00 | fail |
| rarefy to min depth + Wilcoxon | 71 | 8 | 63 | 0.11 | 1.00 | fail |
| log-CPM + OLS(sex) | 63 | 8 | 55 | 0.13 | 1.00 | fail |
| NB GLM, log(library size) offset | 100 | 8 | 92 | 0.08 | 1.00 | fail |

Threshold rationale: every naive method has precision <= 0.13 and recall 1.0 (it finds all 8 plus dozens of
false hits), so any precision threshold in 0.2-0.7 discriminates; 0.40 was chosen at the low end so that legitimate
but sloppier compositional pipelines pass (pseudocount 1 with no filtering reaches 0.50 because pseudocounts bias
sparse taxa in CLR space). 0.75 recall tolerates missing 2 of the 3 modest decliners.

Seed sweep (20 seeds, same generator): compositional methods call 7-15 taxa with 7-8 of 8 truth (precision >= 0.67 in
all seeds for the reference `correct`, worst case 12 called / 8 TP); Wilcoxon-on-proportions calls 36-90 (precision <= 0.22).
The frozen seed is the plain constant, not the best of the sweep.

## Findings and blind spots (honest)
- **Pure CLR passes.** It is compositional; the CLR centre (geometric mean) shifts only ~0.04 log units here because
  only 8 of ~185 taxa change. CLR would be biased if a large fraction of taxa changed asymmetrically. Not tested here.
- **DESeq2 / edgeR-TMM style normalisation passes** (median-of-ratios / trimmed means assume most taxa unchanged, i.e.
  they are reference-frame centring). Only *total-count* normalisation fails. The NB-offset trap uses library size on
  purpose. An agent that runs DESeq2 default will pass; that is correct behaviour, not a leak.
- Pseudocount + no prevalence filter degrades CLR (FPs among very sparse taxa and a few abundant ones) but still passes.
  Pseudocount 1 with more sparse taxa or weaker thresholds could fail an otherwise sound pipeline. Grader margin is 0.10.
- Sex adjustment is nearly irrelevant to the outcome; no compositional method calls a sex-effect taxon (Wilcoxon on proportions calls 4, but it fails anyway). So the
  "adjust for sex if appropriate" clause is inert (does not test covariate handling).
- ALDEx2 / ANCOMBC / LinDA / edgeR / DESeq2 are NOT installed in the sandbox image (only MASS, lme4 among those checked).
  A competent agent must implement CLR / centring by hand in Python or R, which is easy.
- A trivially cheap gaming route exists: answer only the top ~8 taxa by |CLR difference| (no testing). Passes (7/8 TP on frozen data), and is
  essentially the same answer as the correct one; acceptable.

## Ambiguity / risks
- The prompt asks "which taxa differ in abundance" and says nothing about absolute vs relative. This is the intended
  trap, but a solver could argue the relative-abundance reading is literal. Data provide no absolute-load measurement, so
  the only defensible target is the underlying abundance change, which is the field convention for DA tools.
- The 5% BH is on the analyst's own test set (prevalence filtering changes the denominator); the grader is
  precision/recall on the list, so this does not matter.
- All-zero taxa yield NaN in some tests (Wilcoxon, OLS); agents may need to drop them. Mild friction, deliberate.
- Answer ids are case-insensitive in the grader.

## Reproduce
`python generate.py && python solvers.py` (~10 s). Generation is byte-deterministic (LF newlines, fixed seed).
