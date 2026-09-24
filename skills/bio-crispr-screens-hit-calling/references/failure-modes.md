# Reconciliation Failure Modes

Read when two or more methods disagree, or a consensus list comes back empty or inflated.

## Failure Modes

### MAGeCK and BAGEL2 disagree by 200+ hits at FDR 0.05

**Trigger:** Heavy-selection screen (>40% guides change) or cancer-line CN bias.
**Mechanism:** MAGeCK median normalization breaks; BAGEL2 is robust due to reference-set anchoring.
**Symptom:** MAGeCK hit list inflated; BAGEL2 list closer to expected size.
**Fix:** Run MAGeCK with `--norm-method control`; apply CN correction; trust BAGEL2 for essentiality.

### Chronos and MAGeCK disagree at the top 10 in a cancer line

**Trigger:** Top hits are at amplified loci.
**Mechanism:** Chronos models CN bias; MAGeCK does not.
**Symptom:** ERBB2 in HER2+, MYC in MYC-amplified, etc. are top hits in MAGeCK but not Chronos.
**Fix:** Apply [[copy-number-correction]] before MAGeCK or switch to Chronos.

### drugZ and MAGeCK disagree on small-effect drug-modifier screen

**Trigger:** Effect size is small; MAGeCK rank-based test is less sensitive than drugZ bidirectional Z.
**Mechanism:** drugZ specifically optimized for small effects in drug screens (Colic et al. 2019); MAGeCK RRA loses sensitivity at small effects.
**Symptom:** At matched FDR, drugZ calls small-effect chemogenomic interactions (e.g. DDR genes) that MAGeCK RRA misses, with stronger expected-pathway enrichment.
**Fix:** Use drugZ as primary for chemogenomic; MAGeCK as confirmatory. See [[drugz-chemogenomic]].

### JACKS down-weights efficiency, MAGeCK doesn't, disagreement

**Trigger:** A gene has one or two strong sgRNAs and 2-3 weak ones; MAGeCK averages them, JACKS down-weights the weak.
**Mechanism:** JACKS variational Bayes correctly identifies low-efficacy guides; MAGeCK aggregates without this prior.
**Symptom:** Gene is JACKS hit but not MAGeCK.
**Fix:** Inspect per-sgRNA LFC; if strong guides are consistent, JACKS is correct. Validate gene orthogonally.

### Consensus across 3 methods is empty (no hits)

**Trigger:** Either no real biology, each method hits a different failure mode, or -- easy to
overlook -- the merged files are not from the same experimental comparison.
**Mechanism:** Screen quality is low and signal-to-noise across all methods is poor; OR each
input file is individually valid but answers a different question (e.g. a real essentiality
MAGeCK+BAGEL2 pair merged against a drugZ table from an unrelated drug-vs-vehicle screen run on
the same library). `consensus_hits()` merges without error either way -- verified on real data:
merging matched MAGeCK/BAGEL2 essentiality output with a real but unrelated drugZ drug-screen
table produced a 0-gene Tier-1 consensus even though the underlying essentiality screen has
100% precision against CEGv2/NEGv1 (see Input 1's numbers).
**Symptom:** Tier 1 consensus list is empty.
**Fix:** Check each input file's own QC first (screen-qc PR-AUC, precision/recall against
CEGv2/NEGv1) -- a high-precision screen with an empty 3-method consensus points to a mismatched
file, not a bad screen. `consensus_hits()`'s `_check_comparable()` warning (in `scripts/consensus_hits.py`) also fires
when a pair's hit-set overlap is no better than chance, the statistical signature of a
mismatched comparison. Only re-audit QC once a same-comparison mismatch has been ruled out.

