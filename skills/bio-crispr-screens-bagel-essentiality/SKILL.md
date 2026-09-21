---
name: bio-crispr-screens-bagel-essentiality
description: Identifies essential genes from CRISPR-Cas9 fitness screens using BAGEL2 (Kim & Hart 2021 Genome Med), a Bayesian classifier scoring per-gene Bayes Factors via log-likelihood ratios over per-sgRNA fold changes, calibrated against CEGv2 core-essentials (Hart 2017 G3, ~684 genes) and NEGv1 non-essentials (Hart 2014, ~927 genes). Covers the fc + bf + pr workflow, the linear-extrapolation improvement over BAGEL1 truncation, multi-target off-target correction, tumor-suppressor sensitivity (BAGEL2 detects enrichment), and BF calibration (BF >6 ≈ 90% posterior per Hart 2017; ~5% FDR by BAGEL convention). Use when classifying essential vs non-essential genes, calibrating BAGEL2 thresholds against PR curves, identifying tumor suppressors alongside essentials, comparing BAGEL2 hits to MAGeCK / drugZ, or generating publication-quality essentiality calls.
tool_type: cli
primary_tool: BAGEL2
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: BAGEL2 2.0 (hart-lab/bagel, build 115), pandas 2.2+, numpy 1.26+, scipy 1.12+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `BAGEL.py fc --help`; `BAGEL.py bf --help`; `BAGEL.py pr --help`
- Python: BAGEL2 is distributed via `git clone` (no canonical PyPI release); confirm `BAGEL.py version` after checkout.

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## BAGEL2 Essentiality Analysis

**"Identify essential genes from my CRISPR fitness screen using BAGEL2"** -> Compute per-sgRNA fold changes from counts, derive per-gene log-likelihood ratios against reference essential and non-essential gene sets, sum to Bayes Factor, and apply BF threshold calibrated by precision-recall against the reference.

- CLI: `BAGEL.py fc` to compute fold changes
- CLI: `BAGEL.py bf` to compute Bayes Factors
- CLI: `BAGEL.py pr` for precision-recall curves
- Reference sets: CEGv2 (essentials) and NEGv1 (non-essentials); both at https://github.com/hart-lab/bagel

## The BAGEL2 Bayesian Framework (under the hood)

**Why this matters for postdoc-level use:** BAGEL2 uses a Bayes-factor classifier trained on known essential and non-essential genes. The chain:

1. For each sgRNA, compute log-fold-change (LFC) treatment vs control.
2. For each gene, look up per-sgRNA LFCs.
3. For each sgRNA, compute the log-likelihood ratio: `log( P(LFC | gene is essential) / P(LFC | gene is non-essential) )`. The numerator and denominator are KDEs (kernel density estimates) of LFC distributions from CEGv2 and NEGv1 reference sgRNAs.
4. Sum per-gene log-likelihood ratios across all sgRNAs targeting the gene -> per-gene Bayes Factor.
5. Resampling for the confidence interval (default: 10-fold cross-validation; `-b` switches to bootstrapping with `-NB`, default 1000); BF >6 corresponds to ~90% posterior probability (Hart 2017 G3); ~5% FDR by BAGEL convention.
6. **Resampling is randomized and BAGEL.py seeds it from the clock by default -- always pass `-s <fixed-int>`.** See "Reproducibility: Fixing the Random Seed" below; this is required for any result described as publication- or clinical-grade.

**Critical BAGEL2 improvements over BAGEL1:**

- **Linear extrapolation**: BAGEL1 truncated the LLR at the edges of its KDE; BAGEL2 fits a linear regression in the stable region and extrapolates, giving wider dynamic range. This recovers tumor suppressors (highly positive LFC) that BAGEL1 missed.
- **Multi-target correction**: For sgRNAs targeting multiple genomic loci (off-targets), BAGEL2 discards them and regresses out their BF contribution, but only when `-m/--filter-multi-target` is given together with `--align-info`. The original BAGEL counted off-target hits as essentiality signal.
- **Tumor suppressor sensitivity**: BAGEL2 correctly identifies positive selection (enrichment) genes -- not possible in BAGEL1.

## Calibration to CEGv2 / NEGv1

**Why these reference sets matter:** BAGEL2's discriminative power depends on KDEs of LFCs from known essential vs known non-essential genes. CEGv2 (Hart 2017) is 684 core essential genes shared across cell lines; NEGv1 (Hart 2014) is 927 non-essential genes verified across multiple screens. These act as positive and negative controls within every screen.

Reference set integrity:
- CEGv2: pan-cancer essentials -- common dropouts across most cancer cell lines
- NEGv1: confidently non-essential -- genes without expression or genes with verified neutral status

**Critical pitfall:** Using a custom essentiality reference (e.g., a single-cell-line CRISPR screen) instead of CEGv2 biases the BAGEL2 model toward that line's specific biology. Always use the standardized references unless there is a specific reason for custom training.

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
# -s: fixed seed -- required for reproducible calls, see "Reproducibility" below
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
triage, you must opt into bootstrap -- `-b -NB 1000` takes roughly 20-25 minutes at
genome scale (18k+ genes); the CV default runs in well under a minute but has no STD
column.

**Interpretation rule:** BF >6 corresponds to ~90% posterior probability of essentiality against CEGv2 (Hart 2017; FDR ≤3% in that calibration, with ~5% a looser BAGEL convention); higher BF = stronger evidence the gene is essential. BAGEL2 also reports negative BFs which can indicate tumor suppressors (positive selection) -- see "Interpret BAGEL2 Results" below for when that call is actually defensible.

## Reproducibility: Fixing the Random Seed

**BAGEL.py's resampling (both cross-validation and bootstrap) is seeded from the clock
by default** (`-s` defaults to `int(time.time()*100000 % 100000)`, undocumented in
SKILL.md, usage-guide.md or `run_bagel2.sh` before this fix) -- every unseeded run
produces different Bayes Factors. Verified on real HAP1 TKOv3 data (build 115): two
identical, unseeded `bf` invocations differed by up to 26.7 BF and **flipped 33 of
18,053 genes across the BF>6 essential/non-essential threshold**; re-running the
comparison with `sort`+`diff`, 18,026/18,053 rows differed between the two runs.

**Fix:** always pass `-s <fixed-integer>` and record the value used for any
publication- or clinical-grade call.

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
| 6 | Standard call (~90% posterior, Hart 2017) |
| 12 | High-confidence |
| 30 | Ultra-stringent; near-certain essentials |

**Pick threshold based on application:** For exploratory hit calling, BF >0 with low precision is acceptable; for clinical-grade essentiality calls, BF >12 or higher.

## Interpret BAGEL2 Results

**Goal:** Stratify genes into essential, non-essential, and tumor-suppressor categories.

**Approach:** Apply BF threshold to classify essentials; only classify negative-BF genes as tumor suppressors when the screen design actually expects enrichment (see Failure Modes below) -- verified on real HAP1 TKOv3 data (a dropout screen) that skipping this guard flags **86.5% of the genome** as tumor-suppressor, with three assay-control pseudo-genes (`LacZ`, `luciferase`, `EGFP`) as the top hits.

```python
import pandas as pd
import warnings

# Assay-control pseudo-genes spiked into CRISPR libraries (never real biology) --
# verified present and dominating the naive tumor-suppressor call on real HAP1 TKOv3
# output; extend this set to match your library's own controls.
ASSAY_CONTROLS = {'LacZ', 'luciferase', 'EGFP'}

def interpret_bagel(bf_path, bf_essential=6, bf_tumor_suppressor=-6,
                     screen_type='dropout', control_genes=ASSAY_CONTROLS,
                     tumor_suppressor_frac_warn=0.05):
    '''Classify genes from BAGEL2 BF output.

    screen_type: 'dropout' (default) only calls `essential`. Tumor-suppressor calls
    require screen_type='enrichment' or 'both' -- per the Failure Modes section below,
    a pure dropout screen's negative-BF genes are noise, not tumor suppressors.
    '''
    df = pd.read_csv(bf_path, sep='\t')
    df = df[~df['GENE'].isin(control_genes)].copy()   # drop assay-control pseudo-genes
    df['call'] = 'neutral'
    df.loc[df['BF'] > bf_essential, 'call'] = 'essential'
    if screen_type in ('enrichment', 'both'):
        df.loc[df['BF'] < bf_tumor_suppressor, 'call'] = 'tumor_suppressor'
        frac = (df['call'] == 'tumor_suppressor').mean()
        if frac > tumor_suppressor_frac_warn:
            warnings.warn(
                f"{frac:.1%} of genes flagged tumor_suppressor -- implausibly high; "
                "this usually means a dropout-only screen is being scored for "
                "enrichment. Re-check screen_type and BF<-6 calls against literature "
                "before reporting."
            )
    return df.sort_values('BF', ascending=False)
```

Verified on real HAP1 TKOv3 data (a T0-vs-T18 dropout screen): the default
(`screen_type='dropout'`) returns 0 tumor-suppressor calls and excludes the 3 control
pseudo-genes; explicit `screen_type='enrichment'` still surfaces the true tumor
suppressors TSC1/TSC2 as the two most-negative-BF entries (now that the control genes
that previously masked them are excluded) but raises the 86.6%-flagged warning so the
caller doesn't report it uncritically.

**Tumor suppressor identification:** Genes with significantly negative BF (e.g., <-6) in a screen actually designed to detect enrichment (drug-resistance, GoF) indicate fitness advantage from their loss. This is biologically distinct from "non-essential" and may indicate tumor-suppressor function. BAGEL1 could not detect this; BAGEL2's linear extrapolation enables it. **Do not call tumor suppressors from a pure dropout screen** -- see Failure Modes.

## Bayesian Reasoning Per Sgrna

**Why this matters:** BAGEL2 computes per-sgRNA contributions; a gene with 4 sgRNAs each contributing +5 to BF gets +20 total. A gene with 3 sgRNAs contributing +5 and 1 sgRNA contributing -3 (off-target or low-efficacy) gets +12 net.

**Approach:** Add `-r/--sgrna-bayes-factors` to the `bf` command -- this flag is not
mentioned elsewhere in this Skill but is required to get per-sgRNA output; found via
`BAGEL.py bf --help`, not documented upstream either.

```bash
BAGEL.py bf \
    -i foldchange.foldchange \
    -o bayes_factor_sgrna.txt \
    -e CEGv2.txt -n NEGv1.txt \
    -c Sample1,Sample2,Sample3 \
    -s 42 -r                               # -r: per-sgRNA BF contributions
# Output columns: RNA  GENE  <sample columns>  BF  (one row per sgRNA)
```

Verified on real HAP1 TKOv3 data: for RPS3 (the example gene below), the 4 per-sgRNA
BF values summed to 78.9 against a gene-level BF of 80.8 (within 2.4%), confirming the
additive-LLR model this section describes.

**Critical:** When per-sgRNA contributions are very heterogeneous (one sgRNA dominates BF), the gene is "guide-of-one"; verify with JACKS efficiency analysis or apply the second-best-sgRNA rule from [[hit-calling]].

## Comparing BAGEL2, MAGeCK, drugZ

| Property | BAGEL2 | MAGeCK | drugZ |
|----------|--------|--------|-------|
| Statistical framework | Bayes factor with reference sets | NB GLM | Bidirectional Z-score |
| Calibrated against | CEGv2 / NEGv1 | Internal null | Vehicle distribution |
| Tumor suppressor detection | YES | Limited (RRA positive-selection score) | YES |
| Best for | Essentiality classification | General hit calling | Chemogenomic drug screens |
| Output | Bayes factor + CI | FDR + LFC | Z-score + FDR per direction |
| Hit threshold | BF >6 | FDR <0.05 | FDR <0.05 |
| Library calibration | Indirect (reference set) | None | None |

**Reconciliation:** BF >6 ≈ 90% posterior probability (Hart 2017 G3); commonly treated as roughly MAGeCK FDR 0.05 by convention. BAGEL2 hits absent from MAGeCK suggest weak signal that BAGEL2's reference anchoring detects but MAGeCK's null-based test misses; verify by inspecting per-sgRNA contributions.

## Failure Modes

### BAGEL2 returns no hits, crashes, or silently corrupts output on a bad reference set

**Trigger:** Wrong reference gene set file -- species mismatch, thin/unrepresentative reference, or `-e`/`-n` swapped.
**Mechanism:** The actual failure mode depends on *how* the reference is wrong, verified directly against three concrete misconfigurations (not one uniform symptom):
- **Species mismatch** (e.g., a mouse CEG file against a human screen): hard crash, `ValueError: dataset input should have multiple elements` (`scipy.stats.gaussian_kde`) -- not a silent "no hits" run.
- **`-e`/`-n` arguments swapped** (an easy transcription mistake, not covered anywhere else in this Skill): `BAGEL.py bf` exits 0 and writes a complete-looking `bayes_factor.txt`, but **every gene's `BF` is the literal string `nan`** -- verified across all 18,053 genes on real data, no warning surfaced to the user.
- **Genuinely thin/unrepresentative reference** (too few genes, weak KDE separation): the "median BF near zero, no genes >6" symptom below does apply here.
**Symptom:** Depends on failure type -- hard crash (species mismatch), all-`nan` BF column (swapped `-e`/`-n`), or median BF near zero with no genes >6 (thin/unrepresentative reference).
**Fix:** Re-download CEGv2 / NEGv1 from https://github.com/hart-lab/bagel. Verify gene symbols match the screen's annotation. Before trusting a "no hits" reading, check `bayes_factor.txt` for an all-`nan` `BF` column -- that means `-e`/`-n` were swapped, not that the screen is flat.

### BAGEL2 calls negative-LFC genes "tumor suppressors"

**Trigger:** Heavy dropout screen where many genes drop out; the dropout signal is captured as positive BF but the *enriched* genes (negative BF) are noise.
**Mechanism:** BAGEL2's symmetric distribution treats deeply enriched genes as significant; in a dropout-only screen, the enrichment signal is purely noise.
**Symptom:** Many genes with negative BF; these don't validate as tumor suppressors.
**Fix:** Restrict tumor-suppressor calling to screens specifically expecting enrichment (e.g., drug-resistance, GoF screens); for dropout screens, only interpret positive BF.

### Bootstrap CI is wide; BF estimates unstable

**Trigger:** Per-gene number of sgRNAs too low (e.g., <4 in some libraries).
**Mechanism:** Bootstrap of LLR over very few sgRNAs creates wide CI.
**Symptom:** STD column larger than BF; many genes have CI spanning zero.
**Fix:** Use a library with at least 4-6 sgRNAs/gene; or switch to bootstrapping (`-b -NB 5000`); or filter out genes with <3 sgRNAs.

### Low BF for known essential despite high LFC

**Trigger:** One sgRNA per gene is contributing very low LLR (off-target or low-efficacy).
**Mechanism:** BAGEL2 sums LLR; one weak guide drags total down.
**Symptom:** Known essential like RPS3 has BF <6 despite 3 of 4 guides showing -5 LFC.
**Fix:** Inspect per-sgRNA LLR; identify the dragging guide; verify whether to exclude or to use JACKS for efficacy-aware analysis.

### Non-cancer cell-line screen with custom essentials

**Trigger:** Human embryonic kidney HEK293T or iPSC-derived neurons where standard essentials may not be essential.
**Mechanism:** CEGv2 is calibrated for cancer cell lines; some essentials in tumor cells are not essential in iPSC.
**Symptom:** PR curve against CEGv2 shows poor separation; many CEGv2 essentials don't drop out.
**Fix:** Use a cell-type-specific essentialome derived for the relevant lineage; or use MAGeCK / Chronos which doesn't depend on reference sets.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Standard essentiality call | BF >6 | Hart 2017: BF>=6 ~ 90% posterior |
| Stricter essentiality call | BF >12 | BAGEL convention; regenerate precision/recall per screen with `BAGEL.py pr` |
| Ultra-stringent call | BF >30 | BAGEL convention |
| BF for tumor-suppressor candidate | <-6 | Empirical; verify with orthogonal screen |
| Resampling | 10-fold cross-validation (default); `-b -NB 1000` to bootstrap | BAGEL2 default |
| Seed for reproducible/publication-grade calls | `-s <fixed-int>`, always | Unseeded runs flip up to 33/18,053 gene calls at BF>6 (verified on HAP1 TKOv3, build 115) |
| Min reads per sgRNA in control | 30 | Convention; the BAGEL2 default is 0 |
| Min sgRNAs per gene for stable BF | 4-6 | Wider with library convention |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| No hits despite essentials present | Wrong reference set | Re-verify CEGv2 / NEGv1 files |
| Every gene's `BF` is `nan` | `-e`/`-n` arguments swapped | Check argument order against the Output columns example; re-run |
| Different BF/hit list on a rerun with identical inputs | No `-s` seed given (BAGEL.py defaults to a clock-based seed) | Always pass `-s <fixed-int>`, see Reproducibility section |
| Wide resampling CI | Too few sgRNAs/gene | Increase library coverage; bootstrap with more iterations |
| Negative BF for known essentials | Confounding factor (e.g., CN amplification) | Pre-correct with CRISPRcleanR / Chronos |
| Tumor suppressor calls don't validate | Pure dropout screen; enrichment is noise | Only call `interpret_bagel(..., screen_type='enrichment')` on screens actually expecting enrichment |
| Per-sgRNA LLR dominated by one guide | Outlier or off-target | Apply second-best-sgRNA rule |

## References

- Kim E & Hart T. 2021. *Genome Medicine* 13:2. BAGEL2 algorithm and improvements.
- Hart T & Moffat J. 2016. *BMC Bioinformatics* 17:164. BAGEL Bayes factor framework.
- Hart T et al. 2017. *G3* 7:2719. CEGv2 core-essential reference set; BF posterior calibration.
- Hart T et al. 2014. *Mol Syst Biol* 10:733. Gold-standard essential and non-essential reference sets; source of NEGv1.
- Pacini C et al. 2021. *Nat Commun* 12:1661. Integrated cross-study dependencies; reference essentiality benchmarks.

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK RRA/MLE alternative
- crispr-screens/jacks-analysis - JACKS for per-guide efficacy
- crispr-screens/drugz-chemogenomic - drugZ for drug screens
- crispr-screens/hit-calling - Cross-method decision tree
- crispr-screens/screen-qc - Pre-BAGEL QC including CEGv2 PR-AUC
- crispr-screens/library-design - 4-6 sgRNAs/gene library standard
- crispr-screens/copy-number-correction - Pre-correction for cancer-line screens
- pathway-analysis/go-enrichment - Downstream functional analysis
