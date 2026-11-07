# i6_proteomics_de (not agent-visible)

## Pitfall and source
Label-free proteomics missing values are largely missing-not-at-random (low-abundance proteins fall below the detection limit), and per-sample loading differs. Naive min/zero imputation plus a t-test distorts variance and means; un-normalised data are biased by loading (Lazar et al. 2016, J Proteome Res 15:1116; standard practice). Designed from that pitfall and general practice, not from any Skill.

## World (generate.py, seed 1, deterministic; verified byte-identical on rerun)
- 1,480 real protein groups + 12 `REV__` + 8 `CON__` rows (flagged `+`), 12 samples (6 treated / 6 control), 2 batches each 3+3. Raw linear intensities, empty cells = missing. ~147 KB.
- log2 y = mu (N(23.3, 2.4)) + per-protein batch-2 shift (sd 0.6) + condition FC + per-sample loading + noise. Noise sd is larger at low abundance with an inverse-chi2 spread (limma-like, d0=8).
- Loading: uniform over 1 log2 (2x), rejection-sampled so treated-minus-control mean offset is 0.30-0.45 (realised -0.334).
- Missingness: p(drop) = expit(-(y - 21.4)/0.55), plus 2% MCAR. Overall 26.5% missing.
- 150 designed DE proteins, |log2FC| U(0.9, 2.0), 50/50 up/down, random positions (low-abundance included).
- 4 contaminant rows carry a large condition difference; 4 of the 14 eligible flagged rows are called significant by the correct method if they are not excluded (so the exclusion instruction matters slightly).
- Truth (`truth.json`, grader `canonical_markers`) = designed-DE proteins with >= 3 valid values in each condition (the eligibility rule in the task text): 107 proteins (64 up, 43 down). 43 of the 150 are lost to eligibility, mostly low-abundance proteins knocked out in treated.

## Grader
`marker_gene_precision_recall`, answer field `differentially_abundant_proteins`, precision >= 0.85, recall >= 0.80 (no `expected_count`). Case-insensitive matching; IDs are unique upper-case.
Calibration: sound methods land at precision 0.93-0.96, recall 0.88-0.94. Named traps have recall 0.43-0.64 and/or precision 0.65-0.83. Thresholds sit in the gap; recall is the discriminator (margin >= 0.08 for the sound OLS variant, >= 0.16 for named traps).

## Results (solvers.py, seed 1, truth = 107)
| solver | n | TP | prec | recall | verdict |
| --- | --- | --- | --- | --- | --- |
| correct: median-norm, eligibility, moderated t (own limma fitFDist) + batch, BH | 108 | 100 | 0.926 | 0.935 | pass |
| correct: same, OLS t + batch | 98 | 94 | 0.959 | 0.879 | pass |
| trap: no norm + per-protein min imputation + t | 82 | 68 | 0.829 | 0.636 | fail |
| trap variant: global-min imputation | 51 | 46 | 0.902 | 0.430 | fail |
| trap: zero imputation (log2(x+1), normalised) | 69 | 68 | 0.986 | 0.636 | fail |
| trap variant: zero imputation, no norm | 50 | 45 | 0.900 | 0.421 | fail |
| trap: Welch, complete cases, raw p<0.05 | 95 | 62 | 0.653 | 0.579 | fail |
| trap: normalise, ignore batch (moderated t) | 95 | 92 | 0.968 | 0.860 | PASS = BLIND SPOT |
| info: Perseus-style downshift imputation + moderated t + batch | 86 | 84 | 0.977 | 0.785 | fail (hair below 0.80) |
| info: norm + per-protein min imp + Student t | 86 | 83 | 0.965 | 0.776 | fail (hair) |
| info: norm + Student t on observed, no batch | 84 | 81 | 0.964 | 0.757 | fail (hair) |

## Findings that differ from the brief's expectation
- The imputation traps fail mainly on **recall (lost power)**, not by manufacturing false positives. With eligibility enforced and n=6, min/zero fill inflates within-group variance; precision stays 0.83-0.99. Only "no normalisation" adds false positives (loading imbalance). Un-normalised + min fill: precision 0.83, but recall is the decisive failure.
- Median-centring under MNAR is slightly biased (samples with more dropouts have an upward-shifted observed median): per-sample error up to 0.2 log2, treated-minus-control 0.06. Harmless at loading imbalance 0.33; at imbalance ~0.5 even the "correct" median-centred analysis reached precision 0.6 in a sweep. That is why the imbalance range was capped at 0.45. A complete-case median would be a better normaliser.

## Known blind spots
- **Ignoring batch passes** (moderated t, 0.968/0.860). The design is balanced, so batch only inflates residual variance (costs power, no bias). Batch sd 0.9 still only dropped recall to 0.71 (checked in a sweep), so it was left at "modest" 0.6. Student t without moderation and without batch fails narrowly (0.757), i.e. the outcome for "ignore batch" depends on the test used, not on batch per se.
- Grader cannot tell why a list is right; a method with the right power passes even without imputation-aware handling (e.g. observed-values-only OLS passes; that is legitimate and intended).
- Perseus-style imputation (also legitimate) fails narrowly on recall (0.785); other imputation seeds varied 0.70-0.83 across data seeds. Do not read it as a trap result.

## Robustness (seeds 1-4 tried for the world; final = seed 1, not selected on outcome)
Moderated correct passes on all four (recall 0.94-0.99, precision 0.87-0.96). OLS+batch failed precision on seed 2 (0.835, loading imbalance -0.44). Named traps failed on all four seeds. The 0.80 recall threshold is therefore not a hard margin for OLS variants on other draws.

## Ways the task could be gamed / ambiguous
- Answer-size gaming: the truth list is 107 of 1,480; listing every eligible protein fails precision.
- Eligibility: reporting a protein with <3 valid values in a condition is counted as a false positive (truth cannot contain it). "Valid" = non-empty cell.
- ID format: use `Protein IDs` (one ID per row, no semicolon groups). Gene names are also unique but the task asks for IDs.
- Empty cells vs zero: the file uses empty cells; some agents may read them as NaN and others as blanks, both fine with pandas/readr.
- Direction/threshold: no fold-change cut-off is imposed, only 5% BH FDR; a solver applying an extra |log2FC| filter would lose recall and could fail.
