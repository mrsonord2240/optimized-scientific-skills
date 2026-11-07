## rMATS-turbo Differential Analysis

**Goal:** Detect statistically significant differential splicing between two groups from BAMs.

**Approach:** Run rMATS-turbo without `--statoff`, then filter by FDR + ΔPSI + per-replicate coverage.

```bash
# -t paired = paired-end reads (use -t single for single-end); --readLength = the read length
rmats.py \
    --b1 condition1_bams.txt \
    --b2 condition2_bams.txt \
    --gtf annotation.gtf \
    -t paired \
    --readLength 150 \
    --variable-read-length \
    --libType fr-firststrand \
    --nthread 8 \
    --od rmats_output \
    --tmp rmats_tmp \
    --novelSS \
    --cstat 0.05
```

`--cstat 0.05` tests `|ΔPSI| > 0.05`; raise to 0.10 for stricter discovery. `--novelSS` enables novel-junction discovery (recommended with STAR 2-pass).

**Match the command to the data.** rMATS exits 0 and writes header-only tables when it cannot use the reads:
- `-t paired` on single-end BAMs gave 0 events; use `-t single` for single-end data.
- `--readLength` = the read length. On 75 nt reads, `--readLength 150` without `--variable-read-length` gave 0 events; with `--variable-read-length` it only sets IncFormLen/SkipFormLen. `scripts/diff_splicing_rmats.sh` reads type and length from the first BAM.
- After every run assert that `SE.MATS.JC.txt` has rows and an `FDR` column before filtering.
- Strandedness: `--libType` must match the library (RSeQC `infer_experiment.py`); `fr-firststrand` on unstranded data cut SE events 226 -> 121.

**`--paired-stats`** (paired designs: entry *i* of the comma-separated `--b1` list pairs with entry *i* of `--b2`) needs the R package **PAIRADISE** in the R that rMATS finds. Without it rMATS logs `no package called PAIRADISE`, exits 0 and writes header-only tables with no FDR column, so run the table check above. Install (public, no login): `git clone https://github.com/Xinglab/PAIRADISE`, install R packages `nloptr doParallel foreach iterators`, then `R CMD INSTALL PAIRADISE/pairadise/src/pairadise_model`.

Progress is written to `<od>/tmp/JC_SE/pairadise_status.txt`. On WSL2, PAIRADISE's default socket cluster hung in 3 of 6 runs, including after progress reached 100%; set a wall-clock timeout and inspect the output table rather than trusting the progress display. For idle workers, force forked workers: put `parallel:::setDefaultClusterOptions(type = "FORK")` in a file and run rMATS with `R_PROFILE_USER=<file>`. Checked on 5 vs 5 pairs (120 SE events): 31 s, 22/24 strong planted events found with correct direction, 0 calls on the null comparison.

```python
import pandas as pd
import numpy as np

se = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')

def has_coverage_in_half_or_more_reps(inc, skip, minimum=10):
    """Require total junction coverage in >= half the replicates of one group."""
    inc_values = [int(v) for v in inc.split(',')]
    skip_values = [int(v) for v in skip.split(',')]
    if not inc_values or len(inc_values) != len(skip_values):
        return False
    covered = sum((i + s) >= minimum for i, s in zip(inc_values, skip_values))
    return covered >= (len(inc_values) + 1) // 2

se['coverage_ok'] = se.apply(
    lambda row: has_coverage_in_half_or_more_reps(row['IJC_SAMPLE_1'], row['SJC_SAMPLE_1']) and
                has_coverage_in_half_or_more_reps(row['IJC_SAMPLE_2'], row['SJC_SAMPLE_2']),
    axis=1,
)

significant = se[
    (se['FDR'] < 0.05) &
    (se['IncLevelDifference'].abs() > 0.10) &
    se['coverage_ok']  # >=10 total reads in >= half of each group's replicates
].copy()
```
