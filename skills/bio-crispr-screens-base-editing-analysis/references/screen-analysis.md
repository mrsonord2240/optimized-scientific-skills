# Screen Analysis: Efficiency Filtering, Bystander Attribution, Hit Calling

Read when parsing CRISPResso2 output after a BE screen: filtering sgRNAs by editing efficiency, attributing signal to target vs bystander edits, and aggregating sgRNA scores to variants.

## Editing Efficiency Filtering (Critical Pre-Hit-Calling)

**Goal:** Drop sgRNAs that do not edit efficiently, since unedited reads represent no biological perturbation.

**Approach:** From CRISPResso2 output, compute target-base-conversion percentage per sgRNA; filter library to sgRNAs with >50% target editing in a pilot or co-screened control.

```bash
python scripts/filter_by_editing_efficiency.py <crispresso_outputs_dir> --target-pos 5 --target-base C --threshold 0.5
```

`scripts/filter_by_editing_efficiency.py` drops sgRNAs that edit less than the threshold of reads at the target position. `--target-pos` is the 1-indexed position WITHIN THE QUANTIFICATION WINDOW, in column order (CRISPResso2 does not emit a literal "Position" column). Real CRISPResso2 2.3.4 layout: rows are nucleotide identity (A/C/G/T/N/-); columns are one per window position, header-labeled with the REFERENCE base (headers repeat, so columns are selected positionally); values are FRACTIONS in [0, 1], not 0-100, despite the filename. The script raises a specific `ValueError` on schema drift or an out-of-range position.

**Convention:** Drop sgRNAs below 50% editing for variant-function screens. A common working split is a 30% editing floor for primary screening and a 50% floor for confirmed hits. Below 30%, the screen has insufficient power; above 70%, results approach saturation editing.

## Bystander Edit Attribution

**Why this matters:** When a sgRNA's editing window contains the target base AND a bystander base, the screen scores the combination. To attribute screen signal to the target variant alone, either (a) include sgRNAs that edit only the target (no bystander) -- often impossible -- or (b) deconvolute via parallel measurements.

**Strategies for variant-by-variant attribution:**

1. **Tile multiple sgRNAs with different bystander patterns:** If 5 different sgRNAs all hit the target base but have different bystanders, common signal across them is target-attributable (Hanna 2021 approach).

2. **Use orthogonal chemistry:** Run the same variant scan with prime editor (no bystanders); cross-validate. See [[prime-editing-screens]].

3. **Bystander stratification:** From CRISPResso2 allele table, partition reads by exact edit pattern (target only, target+bystander_1, target+bystander_2, etc.); separately score each pattern's contribution to the phenotype.

4. **Restrict library:** Use only sgRNAs with zero bystanders in the editing window (rare; may exclude most candidate spacers).

```bash
python scripts/deconvolute_bystander.py <CRISPResso_on_x>/Alleles_frequency_table.zip --target-pos 67 --bystander-pos 69
```

`scripts/deconvolute_bystander.py` partitions reads by edit pattern at the target and bystander positions and sums `%Reads` per pattern (the real `Alleles_frequency_table.zip` has no `Reference_pct` column; the per-allele read-fraction column is `%Reads`). Positions are 1-indexed within the aligned allele string.

## Hit Calling for Variant-Function Screens

**Goal:** Score per-variant fitness from a base-editor screen.

**Approach:** Filter library to efficiency-passing sgRNAs (>50% editing), then run MAGeCK MLE or drugZ on the sgRNA-level counts; map each significant sgRNA to its predicted variant + bystander pattern; aggregate to per-variant scores.

```bash
python scripts/aggregate_variant_scores.py <mageck>.sgrna_summary.txt variant_annotation.tsv --out-prefix variant_scores
```

`scripts/aggregate_variant_scores.py` aggregates sgRNA-level scores to per-variant scores. `variant_annotation.tsv` maps each sgRNA to its predicted variants (target + bystanders) and must carry the columns `sgrna`, `target_variant`, `n_bystanders`; MAGeCK's real `sgrna_summary.txt` column is lowercase `sgrna`, so build the annotation with that name. Writes `<prefix>.target_only.tsv` (sgRNAs with no bystanders: mean/std/count of LFC per variant) and `<prefix>.mixed.tsv` (sgRNAs with bystanders).
