## RRA vs MLE Decision Tree

| Experimental design | Recommended | Why |
|---------------------|-------------|-----|
| Two conditions (e.g. drug vs vehicle, treated vs untreated), single cell line, no covariates | `mageck test` (RRA) | RRA is more robust to outlier sgRNAs; faster; default for most published screens |
| Time series (Day 0 -> Day 7 -> Day 14 -> Day 21) | `mageck mle` | RRA cannot model multiple timepoints jointly; MLE estimates per-condition beta scores |
| Multi-cell-line panel (e.g. DepMap-style 5-50 lines) | `mageck mle` with cell-line covariate, or Chronos | MLE handles >2 conditions; Chronos preferred at DepMap scale |
| Paired, two-condition samples (each replicate matched donor/cell prep) | `mageck test --paired` | Use the matching treatment/control order; this is the simpler paired comparison when no additional covariates are needed |
| Paired samples with additional conditions or covariates | `mageck mle` with donor and condition columns | MLE accommodates pairing together with time, dose, cell-line, or batch terms; do not collapse the donor structure |
| Combinatorial (treatment x cell line x time) | `mageck mle` with full factorial design | RRA only handles 1 factor |
| Drug screen (vehicle vs drug, multiple doses) | `mageck mle` with dose covariate OR drugZ (preferred for chemogenomic) | drugZ optimized for chemogenomic; see bio-crispr-screens-drugz-chemogenomic |
| Essentiality (Day 0 -> endpoint, single condition) | `mageck test` for simple dropout; `mageck mle` if multi-cell-line | RRA suffices; or BAGEL2 for Bayesian essentiality; see bio-crispr-screens-bagel-essentiality |
| Multi-batch / multi-screen joint analysis | `mageck mle` with batch covariate, JACKS for guide efficacy, or Chronos | See bio-crispr-screens-batch-correction |

**Fails when:**
- Using RRA for time-series and treating each timepoint as a separate test: false-discovery inflation from un-modeled temporal correlation. Use MLE.
- Using MLE without explicit design matrix in a screen with severe batch effects: beta scores include batch variance. Add batch covariate.
- Using MAGeCK for drug screens without vehicle control: comparing drug vs Day 0 conflates drug effect with proliferation. Compare drug vs vehicle, not Day 0. See bio-crispr-screens-drugz-chemogenomic.

## The RRA Algorithm (under the hood)

**Why this matters for postdoc-level use:** RRA is not "non-parametric"; it tests whether per-gene sgRNA p-value ranks are more clustered toward the extremes than expected under the uniform null. The chain:

1. `mageck count` normalizes raw counts (median by default) and outputs `*.count_normalized.txt`.
2. `mageck test` computes a per-sgRNA p-value under a NB model fitted to per-gene variance (mean-variance trend stabilized via empirical Bayes shrinkage; same family as edgeR/DESeq2 but simpler dispersion estimation).
3. Per-sgRNA p-values are ranked; ranks are normalized to percentile (rho).
4. For each gene with k sgRNAs, the alpha-RRA score is the minimum over `Pr(beta(i, k-i+1) <= rho_i)` for i in 1..k -- the probability of observing the i-th smallest rank in a uniform sample.
5. Gene-level p-value is computed by permuting sgRNA-to-gene assignment to derive an empirical null over alpha-RRA scores.
6. Multiple-testing correction is Benjamini-Hochberg.

**Critical assumption:** RRA assumes most sgRNAs are non-changing (used to estimate the NB dispersion). If >40% of sgRNAs change, median normalization fails and the dispersion estimate is biased. Symptom varies by screen: some heavy-selection screens do call every gene significant, but a real 45%-guides-changing test (see Failure Modes below) instead showed asymmetric recall collapse -- one selection direction masked to 0% recall while the other stayed correct. Both manifestations share the same root cause and the same fix: use control-sgRNA normalization (`--norm-method control`) with non-targeting controls as the reference.

## The MLE Model (under the hood)

**Why this matters:** MLE assumes the per-sgRNA log-count is Negative-Binomial with mean determined by a linear combination of condition betas plus an sgRNA-efficiency term (if supplied). The chain:

1. Define a design matrix where rows are samples and columns are conditions; entries are 1 if the sample belongs to the condition, 0 otherwise. A baseline column (all 1s) is required.
2. Per-gene, the model is `log(count_ij) = baseline_i + sum_c (beta_gc * design_jc) + log(sgrna_efficiency_i) + log(size_factor_j)` where i is sgRNA, j is sample, c is condition.
3. The optimizer finds the per-gene beta_gc maximizing the NB likelihood; per-gene Wald-statistic gives a p-value per condition.
4. Each beta represents the log-fold-change in that condition relative to baseline, accounting for sgRNA efficiency.

**Critical assumption:** Beta scores are interpretable only when the design matrix is correctly specified. A common error is omitting batch as a covariate -- the resulting betas absorb batch variance. The fix is to add batch columns; the betas then estimate biology after batch adjustment.

**Convergence failure:** MLE NaN beta scores indicate optimizer divergence -- usually because a gene has too few sgRNAs with non-zero counts in the relevant conditions. The output includes such genes with NaN; do not interpret them as zero effect.


## Normalization Decision

| `--norm-method` | When to use | Mechanism | Fails when |
|------------------|-------------|-----------|------------|
| `median` (default) | Standard screen, <40% guides change | Scale each sample to a common median of all guides | Heavy selection (>40% guides change) inflates median, biasing scaling |
| `total` | Equal sequencing depth assumed, no outliers | Scale to total reads | Sensitive to PCR jackpots / outlier high-count guides |
| `none` | Already-normalized inputs (DEPRECATED for raw FASTQ counting) | Skips scaling | Almost never appropriate; only for already-normalized inputs |
| `control` | Heavy selection screens; library has ≥500 NTCs | Scale each sample so non-targeting controls have constant median | Requires `--control-sgrna ntcs.txt` listing NTC sgRNAs |

**Diagnostic:** Run with `median` first. If essentialome PR-AUC against CEGv2 is high and Gini is in range, accept. If PR-AUC <0.5 despite reasonable Gini, retry with `control` -- this isolates whether median normalization was masking the signal.

## MAGeCK Test (RRA for Two-Condition)

**Goal:** Identify genes significantly enriched or depleted between treatment and control.

**Approach:** Compute per-sgRNA NB-model p-values, rank, apply alpha-RRA to get per-gene scores, permute for gene-level FDR. Outputs separate columns for positive and negative selection.

```bash
mageck test \
    --count-table screen.count.txt \
    --treatment-id Drug_r1,Drug_r2 \
    --control-id Veh_r1,Veh_r2 \
    --norm-method median \
    --output-prefix drug_vs_veh \
    --gene-lfc-method median \                     # alternative: mean (less robust)
    --sort-criteria pos                            # rank for positive selection (choices: neg, pos; default neg)
# Outputs: drug_vs_veh.gene_summary.txt (gene-level), drug_vs_veh.sgrna_summary.txt
```

**Gene-summary columns:**

| Column | Meaning |
|--------|---------|
| `id` | Gene symbol |
| `num` | Number of sgRNAs in library for this gene |
| `neg|score`, `pos|score` | alpha-RRA score, negative- or positive-selection direction |
| `neg|p-value`, `pos|p-value` | Permutation p-value |
| `neg|fdr`, `pos|fdr` | BH-corrected FDR |
| `neg|rank`, `pos|rank` | Rank by score (1 = top hit in that direction) |
| `neg|lfc`, `pos|lfc` | Median log-fold-change of sgRNAs in that direction |

**Interpretation rule:** A gene is a strong essential if `neg|fdr < 0.05` AND `neg|lfc < -1`. A gene is a strong resistance hit if `pos|fdr < 0.05` AND `pos|lfc > 1`. Genes with `fdr < 0.05` but `|lfc| < 0.5` are statistically significant but biologically weak -- worth flagging for orthogonal validation.

## MAGeCK MLE (Multi-Condition)

**Goal:** Estimate gene effects across complex experimental designs with multiple conditions, time points, batches, or covariates.

**Approach:** Specify a design matrix mapping samples to conditions, then run `mageck mle` which fits a NB GLM with sgRNA-efficiency term and outputs per-gene per-condition beta scores.

```bash
# Design matrix: design.txt (tab-separated)
# Samples must match sample-label in mageck count output
# baseline column must be present and all 1
cat > design.txt <<EOF
Samples	baseline	day7	day14	day21
Day0	1	0	0	0
Day7_r1	1	1	0	0
Day7_r2	1	1	0	0
Day14_r1	1	0	1	0
Day14_r2	1	0	1	0
Day21_r1	1	0	0	1
Day21_r2	1	0	0	1
EOF

mageck mle \
    --count-table screen.count.txt \
    --design-matrix design.txt \
    --output-prefix timecourse_mle \
    --norm-method median \
    --permutation-round 10 \                       # default is 2; raise for any FDR call trusted near 0.05 (see below)
    --sgrna-efficiency efficiency.txt \            # OPTIONAL: from JACKS or library design
    --sgrna-eff-name-column 0 \                    # 0-based; 0 is the default
    --sgrna-eff-score-column 1 \                   # 0-based; 1 is the default
    --max-sgrnapergene-permutation 40              # skip genes with more sgRNAs than this (default 40)
# Outputs: timecourse_mle.gene_summary.txt with beta scores per condition + Wald p-values
```

**Gene-summary columns (MLE):**

| Column | Meaning |
|--------|---------|
| `Gene` | Gene symbol |
| `sgRNA` | Number of sgRNAs |
| `<condition>|beta` | Effect-size estimate (log-fold-change relative to baseline) |
| `<condition>|z` | Wald z-statistic |
| `<condition>|p-value` | Two-sided p-value |
| `<condition>|fdr` | BH-corrected FDR |
| `<condition>|wald-fdr` | Wald-statistic-based FDR (alternative) |

**Interpretation rule:** Beta scores are log2-fold-changes; a beta of -1 in condition day21 means sgRNAs are 2-fold depleted in day-21 relative to baseline. NaN betas indicate convergence failure (typically a gene with too few non-zero counts in that condition); exclude from interpretation.

**Permutation-Round Sensitivity of the primary `fdr` column:** `mageck mle`'s permutation-based `fdr`
column is coarser than it looks -- `--permutation-round` defaults to 2 (`mageck mle --help`: "Suggested
value: 10 (may take longer time). Default 2"). At that default, rerunning the *identical* command can
move near-threshold genes across FDR<0.05: MAGeCK 0.5.9.5 offers no `mle --seed` option to make this
Monte-Carlo component repeatable. Always use an explicitly stated higher round count for a result you
will report, and rerun borderline calls before interpreting them.

One real HAP1 TKOv3 T0-vs-T18 exercise (a 1500-gene/5880-sgRNA random subset, two-condition
`baseline`+`treatment` design, `--norm-method median`) moved from 48 to 67 FDR<0.05 genes between
rounds 2 and 5, with 19 genes changing status, while `wald-fdr` remained 108 across rounds 1, 2, and
5. This is an *illustrative run*, not an expected effect size: a second independently sampled subset
changed by only 3 genes over the same comparison. **Rule of thumb:** don't trust a permutation `fdr`
call near 0.05 at the default round count; use `--permutation-round 10+` (per `--help`'s suggestion),
then check whether the conclusion is stable and cross-check `wald-fdr`. `wald-fdr` is an asymptotic
alternative, not a replacement for validation. If the columns disagree, treat the gene as borderline
and validate orthogonally rather than trusting either one alone.

## Sample MAGeCK Test for Drug Screen with sgRNA Efficiency

```bash
mageck test \
    --count-table screen.count.txt \
    --treatment-id Drug_r1,Drug_r2,Drug_r3 \
    --control-id Veh_r1,Veh_r2,Veh_r3 \
    --norm-method control \                         # NTCs as normalization reference
    --control-sgrna ntcs.txt \                      # one NTC sgRNA per line
    --gene-lfc-method median \
    --variance-estimation-samples Day0_r1,Day0_r2 \  # estimate variance from these samples
    --output-prefix drug_screen_normalized
# --sgrna-efficiency / --sgrna-eff-name-column / --sgrna-eff-score-column are `mageck mle` options,
# not `mageck test` options; pass them to mle if guide-efficacy weighting is needed.
```

**Reading order:** Run `mageck count` -> screen-qc skill for QC -> `mageck test` or `mageck mle` -> `MAGeCKFlute` for visualization -> downstream pathway analysis.


## Failure Modes

### NaN beta scores in MLE output

**Trigger:** A gene has fewer than 2 non-zero-count sgRNAs in the relevant condition.
**Mechanism:** MLE optimizer cannot fit beta when likelihood is flat (insufficient evidence).
**Symptom:** `mle.gene_summary.txt` shows NaN in specific gene-condition cells.
**Fix:** Filter these genes from downstream interpretation; they are not "zero effect" but undetermined. If many genes show NaN, the screen depth is too low or many guides have failed -- audit with screen-qc.

### Heavy-selection screen distorts RRA calls under median normalization

**Trigger:** >40% of sgRNAs change direction (heavy selection screen).
**Mechanism:** Median normalization assumes most guides are non-changing; under heavy selection, median is biased and dispersion estimation breaks.
**Symptom:** Two manifestations of the same root cause, either can occur: (a) thousands of "hits" at FDR <0.05 with a volcano plot that looks like a U with no separation, or (b) one selection direction goes to near-zero recall while the other direction stays correct -- measured directly on a real 200-gene/45%-changing synthetic screen: median normalization gave 100% recall of depleted genes but 0% recall of enriched genes (0/43), while `--norm-method control` recovered both (100%/100%). Don't assume "every gene significant" is the only tell; check recall by direction, not just hit count.
**Fix:** Switch to `--norm-method control` with non-targeting controls; or use BAGEL2 / Chronos for essentiality analysis (`mageck` is not designed for screens with high-fraction true essentiality).

### MLE beta absorbs batch variance

**Trigger:** Multi-batch screen run without explicit batch covariate in design matrix.
**Mechanism:** Without a batch column, beta_condition includes both biological signal and batch difference between conditions.
**Symptom:** Beta scores differ between batches for the same biological condition; PCA shows samples cluster by batch not condition.
**Fix:** Add batch covariate columns to design matrix (e.g., one column per batch indicator); the biological betas are then estimated after batch adjustment. See bio-crispr-screens-batch-correction.

### sgRNA-efficiency injection makes scores worse

**Trigger:** Supplied JACKS efficiency scores from a different cell line / chemistry than the current screen.
**Mechanism:** sgRNA efficacy is cell-line and modality dependent; HEK293T efficiency is not HCT116 efficiency; Cas9 efficiency is not CRISPRi efficiency.
**Symptom:** Some genes that were hits without efficiency become non-significant with efficiency.
**Fix:** Use efficiency derived from JACKS run on the matching cell line / library / chemistry, or use no efficiency (treat all sgRNAs as equal).

### Lack of cell-line covariate in multi-line screen

**Trigger:** Running `mageck mle` on multiple cell lines without a cell-line indicator column.
**Mechanism:** Each cell line has different essentiality profile; combining without indicator pools variance and dilutes per-line signal.
**Symptom:** Per-line known essentials don't show up; results look like a meta-analysis without effect sizes.
**Fix:** Add cell-line indicator columns; or move to Chronos which models cell-line and screen-quality jointly (see bio-crispr-screens-hit-calling).

## Reconciliation: When MAGeCK Disagrees With Other Tools

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MAGeCK FDR significant, BAGEL2 Bayes factor not | Different reference set; BAGEL2 trained on CEGv2/NEGv1 | Trust BAGEL2 for essentiality calling; MAGeCK for general LFC |
| MAGeCK significant, JACKS not | JACKS down-weighted noisy guides; MAGeCK trusts all 4 | Inspect per-sgRNA scores; if one sgRNA dominates, JACKS is right |
| MAGeCK and Chronos agree at top 100; disagree at 100-300 | Chronos accounts for CN and screen quality; MAGeCK does not | Trust Chronos for cancer-line screens; MAGeCK lacks CN correction |
| MAGeCK significant, drugZ not | drugZ uses bidirectional Z; MAGeCK uses RRA | For chemogenomic, trust drugZ (vehicle-anchored) |
| RRA and MLE disagree on same data | RRA more robust to outliers; MLE more sensitive | Higher-ranked hits in both = high confidence; only-one-method = orthogonal validate |

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Hit FDR (gene-level) | <0.05 | MAGeCK paper Li 2014; standard |
| Hit LFC magnitude | >1 (2-fold) | Biologically interpretable effect size; FDR-only hits at low LFC need validation |
| Normalization fraction-changing threshold | <40% guides change for median norm | Median-normalization assumption (most guides unchanged) |
| sgRNAs needed for reliable MLE | ≥3 per gene with non-zero counts in each condition | MAGeCK 0.5+ behavior; below this, NaN |
| Time-course conditions for MLE | ≥3 (otherwise use test) | MLE statistical power |
| PR-AUC of ranked hits against CEGv2 | >0.7 for "passing" essentiality screen | Community convention (CEGv2 from Hart 2017); see bio-crispr-screens-screen-qc |
| RRA permutation passes per gene | 100 default; raise via `--additional-rra-parameters "--permutation N"` | Not exposed directly on `mageck test`; trades runtime |
| MLE `--permutation-round` | Default 2; use >=10 for any reported or near-0.05 FDR call | `mageck mle --help`'s own suggestion; identical default-round reruns can flip boundary calls, and the size of a round-count effect is screen-dependent (see MLE section) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `mageck count` outputs mostly zero counts | Library column order swapped, or wrong `--trim-5` length | Library must be sgRNA id, sequence, gene; try `--trim-5 AUTO` |
| All genes significant, or recall collapses in one selection direction | Median normalization break (heavy selection, >40% guides changing) | `--norm-method control` |
| `FluteRRA`/`FluteMLE` errors on missing output dir, or `undefined columns selected` | `proj=` not passed (subfolder `MAGeCKFlute_<proj>/` not pre-created), or a MAGeCKFlute build ahead of/behind the Version Compatibility note | Pass `proj=` and `dir.create()` the subfolder first; if the column error persists, it is a known MAGeCKFlute-build compatibility gap -- check `packageVersion('MAGeCKFlute')` against this Skill's tested version |
| `FluteMLE`: `Sample name doesn't match` | `ctrlname` or `treatname` was the implicit `baseline` intercept, or is absent from the MLE output | Use two actual emitted condition names (for example `day7` and `day21`); inspect `colnames()` of the MLE result rather than using `baseline` |
| NaN beta in MLE | Gene with insufficient non-zero counts | Exclude from interpretation |
| Two-condition MLE works but ranks differ from RRA | Different test statistic | Both correct; check direction and use the appropriate one |
| Hits include amplified genes | No CN correction | See bio-crispr-screens-copy-number-correction |
| Library not detected in screen.countsummary.txt | Library file format wrong | Check column order is sgRNA id, sequence, gene (no spaces) |
