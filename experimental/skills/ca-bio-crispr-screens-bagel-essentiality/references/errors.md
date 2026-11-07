## Failure Modes

### BAGEL2 returns no hits, crashes, or silently corrupts output on a bad reference set

**Trigger:** Wrong reference gene set file -- species mismatch, thin/unrepresentative reference, or `-e`/`-n` swapped.
**Mechanism:** The actual failure mode depends on *how* the reference is wrong, verified directly against three concrete misconfigurations (not one uniform symptom):
- **Species mismatch** (e.g., a mouse CEG file against a human screen): hard crash, `ValueError: dataset input should have multiple elements` (`scipy.stats.gaussian_kde`) -- not a silent "no hits" run.
- **`-e`/`-n` arguments swapped** (an easy transcription mistake, not covered anywhere else in this Skill): `BAGEL.py bf` exits 0 and writes a complete-looking `bayes_factor.txt`, but **every gene's `BF` is the literal string `nan`** -- verified across all 18,053 genes on real data, no warning surfaced to the user.
- **Genuinely thin/unrepresentative reference** (too few genes, weak KDE separation): the "median BF near zero, no genes >6" symptom below does apply here.
**Fix:** Run the pre-flight and post-run checks above (`scripts/check_bagel_inputs.py`); re-download CEGv2 / NEGv1 from https://github.com/hart-lab/bagel and verify gene symbols match the screen's annotation. An all-`nan` `BF` column means `-e`/`-n` were swapped, not that the screen is flat.

### BAGEL2 calls negative-LFC genes "tumor suppressors"

**Trigger:** Heavy dropout screen where many genes drop out; the dropout signal is captured as positive BF but the *enriched* genes (negative BF) are noise.
**Mechanism:** BAGEL2's symmetric distribution treats deeply enriched genes as significant; in a dropout-only screen, the enrichment signal is purely noise.
**Symptom:** Many genes with negative BF; these don't validate as tumor suppressors.
**Fix:** `scripts/interpret_bagel.py --screen-type dropout` (the default; `references/interpret-calls.md`) makes no tumor-suppressor calls; use `'enrichment'` only for screens expecting positive selection (drug-resistance, GoF); in dropout screens interpret only positive BF.

### Thin per-gene coverage; BF estimates unstable

**Trigger:** Per-gene number of sgRNAs too low (e.g., <4 in some libraries).
**Mechanism:** The LLR is summed over very few sgRNAs and the linear extrapolation amplifies noise.
**Symptom:** Implausibly large BF magnitudes -- verified on a synthetic 3-sgRNA/gene library (200 real CEGv2/NEGv1 genes, 600 sgRNAs, build 115): top BF ~2,400 against a maximum of ~132 on the real 70,754-sgRNA HAP1 TKOv3 library. A wide bootstrap CI (STD larger than |BF|, CI spanning zero) is possible but **not a reliable tell**: 0/200 genes showed it in that run because the planted effects were strong.
**Fix:** Use a library with at least 4-6 sgRNAs/gene; filter out genes with <3 sgRNAs; run `-b` (`-NB 1000` or more) to get the `STD` column for triage, and treat extreme BF magnitudes as suspect even when `STD` looks small.

### Low BF for known essential despite high LFC

**Trigger:** One sgRNA per gene is contributing very low LLR (off-target or low-efficacy).
**Mechanism:** BAGEL2 sums LLR; one weak guide drags total down.
**Symptom:** Known essential like RPS3 has BF <6 despite 3 of 4 guides showing -5 LFC.
**Fix:** Inspect per-sgRNA LLR (`references/per-sgrna-contributions.md`); identify the dragging guide; verify whether to exclude or to use JACKS for efficacy-aware analysis.

### Negative BF for known essentials in a cancer-line screen

**Trigger:** Copy-number amplification confounds the screen (multiple cuts in amplified regions look like dropout, so amplified regions appear "essential" with high BF).
**Fix:** Pre-correct with CRISPRcleanR / Chronos (see copy-number-correction) before BAGEL2.

### Non-cancer cell-line screen with custom essentials

**Trigger:** Human embryonic kidney HEK293T or iPSC-derived neurons where standard essentials may not be essential.
**Mechanism:** CEGv2 is calibrated for cancer cell lines; some essentials in tumor cells are not essential in iPSC.
**Symptom:** PR curve against CEGv2 shows poor separation; many CEGv2 essentials don't drop out.
**Fix:** Use a cell-type-specific essentialome derived for the relevant lineage (the BF >6 calibration from Hart 2017 may not hold outside cancer lines); or use MAGeCK / Chronos which doesn't depend on reference sets.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Essentiality call | BF >6 / >12 / >30 | See the BF ladder under "Precision-Recall Curve" |
| BF for tumor-suppressor candidate | <-6 | Empirical; enrichment screens only; verify with orthogonal screen |
| Resampling | 10-fold cross-validation (default); `-b -NB 1000` to bootstrap | BAGEL2 default |
| Seed | `-s <fixed-int>`, always | See Reproducibility |
| Min reads per sgRNA in control | 30 | Convention; the BAGEL2 default is 0 |
| Min sgRNAs per gene for stable BF | 4-6 | Wider with library convention |
