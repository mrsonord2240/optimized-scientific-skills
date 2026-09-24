# QC report and exclusion-decision template

Fill this in once per dataset, after the raw-signal checks and before normalizing. It records what was
checked, which samples were excluded and why, and whether the downstream result depends on the call.
Copy the file next to the data; do not edit this copy.

Thresholds are not restated here. Each row cites either the `flag` / `investigate` / `possible_swap` /
`status` column a `SKILL.md` function returned, or the row of `SKILL.md` "Quantitative Thresholds" that
applies. A blank cell means the step was not done; it never means clean.

## 1. Dataset

| | |
| --- | --- |
| Study / dataset id | |
| Date, analyst | |
| Acquisition (DDA / DIA / TMT-plexes) and search tool + version | |
| Un-normalized column used for the loading check | (MaxQuant `Intensity <sample>`, DIA-NN `Precursor.Quantity`, TMT raw reporter intensities) |
| Samples per group x batch | (paste the design table; groups of one and fully confounded batches go in section 5) |

## 2. Checks run

One row per check in `SKILL.md`. `not measurable` needs a reason; it is not a pass.

| Check | Function / tool | Result: passed / flagged / not measurable | Why not measurable, or where the flags are logged |
| --- | --- | --- | --- |
| Loading and ID count | `raw_sample_qc` (`loading_rule` says within-group or fallback) | | |
| Contaminant fraction | `contaminant_fraction`, judged against this lab's baseline for the sample type | | |
| Replicate correlation | `replicate_correlation` (summarize `status == 'measured'` only) | | |
| Sample swap | `cross_group_correlation` (`possible_swap`) | | |
| CV | linear-scale CV | | |
| Missingness (MNAR vs MCAR) and completeness | | | |
| PCA / batch | `pca_batch_check` (per-PC `tested` / `not_testable`) | | |
| TMT channel balance | within-plex `fold_vs_plex_median` (`investigate`) | | |
| Level-1 run metrics (DIA-NN) | `diann_level1` (`flag`) | | |
| Queue drift | Levey-Jennings on interspersed QC injections | | |

## 3. Exclusion and decision log

One row per flagged sample or channel, kept or excluded. Stop and ask before excluding (`SKILL.md`, PCA and
batch section): the last column records who approved.

| Sample / channel | Group, batch, plex | Level | Metric (column name) | Value | Threshold (source) | Decision | Reason | Confirmed against | Approved by |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| | | | | | | keep / keep-flagged / exclude / re-run | | sample sheet, re-injection, raw file | |

Decisions:
- **exclude**: a loading or injection failure. It is removed, never rescaled.
- **re-run**: the sample is worth re-injecting or re-preparing before any call.
- **keep-flagged**: kept, and reported as flagged, with its sensitivity row in section 4.
- **keep**: a flag reviewed and judged not a problem, with the reason written.

A swap or relabel candidate (`possible_swap`) is confirmed against the sample sheet before anything is
excluded, because the flag is a candidate, not proof.

Example row (from `examples/qc_analysis.py`, whose seeded run flags `ctrl_3`): sample `ctrl_3`, group
`ctrl`, level 3, metric total raw signal vs the group median (`fold_vs_median` in that script,
`fold_total_vs_group` from `raw_sample_qc`), value 0.30 with 287 of 400 proteins quantified, threshold the loading rule in `SKILL.md` "Inspect Raw Signal and Remove
Contaminants Before Normalizing", decision exclude (loading failure, not rescaled).

## 4. Sensitivity: does the result depend on the call?

Run the downstream check with and without every excluded or keep-flagged sample. Report both, side by side.

| Sample(s) removed | Downstream check | With | Without | Robust? |
| --- | --- | --- | --- | --- |
| | e.g. proteins significant at the chosen FDR; the effect size of the headline protein; PC1 driver; overlap of the top-N list | | | yes / no, and how it changes |

If the conclusion changes with the sample, say so in the summary. Do not report only the version that reads
better.

## 5. Design limits

State what the design cannot show, in words:
- groups with one sample (replicate correlation, CV and the within-group loading rule are not measurable there);
- a batch fully confounded with condition (non-identifiable: report it and stop, do not correct or test);
- fewer than 5 per group (PCA is unstable);
- no un-normalized column available (the loading check cannot run).

## 6. Result

| | |
| --- | --- |
| Samples retained for the differential test | |
| Samples excluded, re-run or kept-flagged | (counts; the rows are in section 3) |
| Normalization applied after exclusion, and how batch enters the design | |
| Overall QC call | passed / passed with flagged samples / not passed |
| Sign-off | |
