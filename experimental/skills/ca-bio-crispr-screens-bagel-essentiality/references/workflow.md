## Compute Per-Sample Fold Changes

**Goal:** Generate per-sgRNA fold-change matrix as input for Bayes-factor calculation.

**Approach:** Take normalized counts, compute log-fold-change vs a control (Day 0 or plasmid baseline) per sgRNA.

```bash
# BAGEL2 installation: distributed via git clone (no canonical PyPI release).
git clone https://github.com/hart-lab/bagel
cd bagel

# Inputs:
# counts.txt: tab-separated with columns: sgRNA, GENE, Sample1, Sample2, ...
# Control column(s): typically Day 0 or plasmid sample(s)
# Treatment column(s): screen endpoint

BAGEL.py fc \
    -i counts.txt \
    -o foldchange \
    -c Plasmid \
    --min-reads 30
# -o is a LABEL for fc, not a file path; writes foldchange.foldchange
# -c: control sample (or Day 0); --min-reads: default is 0, 30 is a common convention
# Output: foldchange.foldchange (per-sgRNA LFCs) and foldchange.normed_readcount
```

**Trap, verified directly: a trailing `# comment` after a line-continuation `\` breaks
the shell.** Comments must go on their own line (as above), not appended after a `\`
that continues to another `-flag` line -- `bash` stops treating the rest of the
pipeline as one command as soon as a `#` follows the escaped space, and every
subsequent `-flag` line then runs as its own (nonexistent) command.

**Pre-flight (run before `bf`):** `BAGEL.py` gives no runtime signal for a species/gene-symbol mismatch (raw traceback) or swapped `-e`/`-n` (exit 0, all-`nan` output). `scripts/check_bagel_inputs.py pre` checks that reference genes are present in the screen and that `-e` genes drop out more than `-n` genes, and exits 1 otherwise:

```bash
python scripts/check_bagel_inputs.py pre foldchange.foldchange CEGv2.txt NEGv1.txt Sample1,Sample2,Sample3
```

## Compute Bayes Factors

**Goal:** Score per-gene essentiality as a Bayes Factor.

**Approach:** Run `BAGEL.py bf` with the fold-change matrix and reference gene sets. Resampling defaults to 10-fold cross-validation; add `-b -NB N` to bootstrap instead.

```bash
BAGEL.py bf \
    -i foldchange.foldchange \
    -o bayes_factor.txt \
    -e CEGv2.txt \
    -n NEGv1.txt \
    -c Sample1,Sample2,Sample3 \
    -s 42 \
    -b -NB 1000
# -e/-n: essentials/non-essentials reference; -c: treatment samples to score
# -s: fixed seed, see "Reproducibility" below
# -b -NB 1000: opt into bootstrapping (default is 10-fold cross-validation)
# Output: bayes_factor.txt - per-gene Bayes Factor + CI
```

**Output columns:**

| Column | Meaning |
|--------|---------|
| `GENE` | Gene symbol |
| `BF` | Per-gene Bayes Factor (log-likelihood ratio summed across sgRNAs) |
| `STD` | Standard deviation across the 10 cross-validation folds (or bootstrap iterations with `-b`) -- **only present when `-b` is given** |
| `NumObs` | Number of sgRNAs contributing -- **only present when `-b` is given** |

Verified directly (build 115): a CV-default run (no `-b`) writes a 2-column `GENE\tBF`
file; only `-b` adds `STD`/`NumObs`. If you need the STD diagnostic for guide-quality
triage, you must opt into bootstrap -- timed directly at genome scale (18,053 genes,
HAP1 TKOv3): ~2.2s/bootstrap iteration, so `-b -NB 1000` takes roughly 35-40 minutes
(longer than an earlier estimate of 20-25 min); the CV default runs in well under a
minute but has no STD column. Runtime is hardware- and load-dependent -- time a short
run (`-NB 50`, a few minutes) first and extrapolate for your own machine.

**Post-run check:** `python scripts/check_bagel_inputs.py post bayes_factor.txt` exits 1 when more than 5% of the `BF` column is `nan` (swapped `-e`/`-n`, species mismatch or wrong `-c`), so an all-`nan` file is never read as "no essentials". Both checks were run on real HAP1 TKOv3 data: they pass on the correct call and fail on swapped `-e`/`-n` (18,053/18,053 `nan`; `-e` mean LFC +0.13 vs `-n` -2.67) and on the mouse CEG file (0/621 genes present).

**Interpretation rule:** BF >6 corresponds to ~90% posterior probability of essentiality against CEGv2 (ladder and FDR mapping under "Precision-Recall Curve"); higher BF = stronger evidence the gene is essential. BAGEL2 also reports negative BFs which can indicate tumor suppressors (positive selection) -- see `references/interpret-calls.md` for when that call is actually defensible.

## Reproducibility: Fixing the Random Seed

**BAGEL.py's resampling (both cross-validation and bootstrap) is seeded from the clock
by default** (`-s` defaults to `int(time.time()*100000 % 100000)`) -- every unseeded run
produces different Bayes Factors. Verified on real HAP1 TKOv3 data (build 115): two
identical, unseeded `bf` invocations differed by up to 26.7 BF and **flipped 33 of
18,053 genes across the BF>6 essential/non-essential threshold**; re-running the
comparison with `sort`+`diff`, 18,026/18,053 rows differed between the two runs.

**Fix:** always pass `-s <fixed-integer>` and record the value used for any
publication-grade call.

```bash
BAGEL.py bf -i foldchange.foldchange -o bayes_factor.txt \
    -e CEGv2.txt -n NEGv1.txt -c Sample1,Sample2,Sample3 \
    -s 42                                   # fixed seed -- record this value
```

Verified directly: two `bf` runs against identical inputs with `-s 42` produced a
byte-identical `bayes_factor.txt` (`diff` clean across all 18,053 genes).

**Trap, verified against build 115: `-s` is declared twice in BAGEL.py's own CLI** --
once for `-s/--use-small-sample` (a flag) and once for `-s/--seed` (an integer).
`BAGEL.py bf --help` prints both and warns `The parameter -s is used more than once`;
passing `-s <int>` binds to `--seed` (the later-declared option wins in Click's
parser) -- confirmed by running `-s 42` and inspecting the resulting fold-assignment
log and byte-identical rerun, not just by reading `--help`. This is BAGEL.py's own
bug, not this Skill's; verify with `--help` before relying on it in an environment
running a different BAGEL2 build.

## Precision-Recall Curve

**Goal:** Empirically select BF threshold for a given precision/recall tradeoff.

**Approach:** Run `BAGEL.py pr` to compute precision and recall at every BF level against CEGv2; pick the BF that gives desired precision.

```bash
BAGEL.py pr \
    -i bayes_factor.txt \
    -o precision_recall.txt \
    -e CEGv2.txt \
    -n NEGv1.txt
# Output: precision_recall.txt - precision/recall at each BF threshold
```

**Practical BF ladder (regenerate precision and recall per screen with `BAGEL.py pr`):**

| BF threshold | Use case |
|--------------|----------|
| 0 | Exploratory; highest recall |
| 6 | Standard call (~90% posterior, Hart 2017; FDR <=3% in that calibration, ~5% a looser BAGEL convention) |
| 12 | High-confidence (BAGEL convention) |
| 30 | Ultra-stringent; near-certain essentials (BAGEL convention) |

**Pick threshold based on application:** For exploratory hit calling, BF >0 with low precision is acceptable; for high-stringency essentiality calls, BF >12 or higher.
