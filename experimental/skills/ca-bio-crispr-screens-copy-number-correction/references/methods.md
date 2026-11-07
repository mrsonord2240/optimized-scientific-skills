## The Copy-Number Artifact (Mechanism)

**Aguirre AJ et al 2016 *Cancer Discov* 6:914** and **Munoz DM et al 2016 *Cancer Discov* 6:900** demonstrated that focal amplification regions in cancer cell lines appear systematically "essential" in CRISPR-Cas9 screens, independent of the gene's actual biology. The mechanism:

1. A focal amplification creates 4-50+ copies of a genomic region.
2. Each sgRNA targeting a gene in that region cuts at all copies simultaneously.
3. Multiple cuts trigger a DNA-damage response and G2 arrest, in both TP53-mutant and TP53-wild-type lines but with larger magnitude in wild-type (Aguirre 2016).
4. Cells arrest in G2 phase; the sgRNA appears depleted because its bearer cells don't proliferate.
5. The depletion is proportional to the number of simultaneous cuts, not the gene's essentiality.

**Consequence:** ERBB2 appears essential in HER2-amplified SK-BR-3. MYC appears essential in MYC-amplified colorectal lines (10+ copies). FGFR1 appears essential in FGFR1-amplified head-and-neck lines. These are all false positives.

**Affects:** All Cas9-KO screens in cancer cell lines. Universal, not conditional. Cannot be remediated by sequencing depth, library size, or replicate count. Requires explicit correction.

The p53-dependence of Cas9-cut toxicity in general was characterized later, by Haapaniemi 2018 and Ihry 2018.

**Bypassed by:**
- CRISPRi (catalytically dead Cas9, no DNA damage) -> no artifact
- CRISPRa (catalytically dead Cas9) -> no artifact
- Base editing (single-strand nick + deaminase) -> reduced artifact
- Prime editing (nick + RT) -> reduced artifact

## Correction Method Decision Tree

| Available data | Recommended method | Why |
|----------------|---------------------|-----|
| Cell-line panel without matched CN profile | CRISPRcleanR | Unsupervised; uses genomic position only |
| Single cell line with matched WGS/SNP-array CN | CRISPRcleanR | Chronos' CN correction (`alternate_CN`) refuses fewer than 3 cell lines, matched CN or not; it says so and points at CRISPRcleanR |
| >=3 cell lines, CN available | Chronos | `alternate_CN`'s minimum; it fits the CN-effect curve across lines |
| DepMap-scale (1000+ cell lines, longitudinal) | Chronos | Population-dynamics + screen quality + CN; DepMap quarterly standard |
| Single cell line, multi-timepoint | Chronos for gene effects, CRISPRcleanR for CN correction | Chronos trains and scores fine on one line, but cannot CN-correct below 3 lines |
| Need to integrate with downstream MAGeCK | CRISPRcleanR (pre-hoc) | Outputs corrected counts for any downstream tool |
| Multiple cell lines + multiple batches | Chronos | Joint modeling of all dimensions |

## CRISPRcleanR (Iorio 2018) - Unsupervised Pre-Hoc

**Goal:** Correct copy-number bias without requiring matched CN profile by detecting position-based systematic enrichment / depletion patterns.

**Approach:** Order sgRNAs by chromosomal coordinate; detect segments where sgRNAs show systematic depletion (or enrichment) inconsistent with single-gene biology; shift these segments toward the global mean. The intuition: focal amplifications create depletion bands extending tens to hundreds of kb; non-amplified essential genes are punctate.

```r
library(CRISPRcleanR)

# Load library annotation (sgRNA -> chromosomal coordinates)
data(KY_Library_v1.0)   # KY library; replace with your library annotation
# OR use ccr.PrepareAnnotations() to make custom

# Load count data with first 2 cols: sgRNA, gene, then sample counts
counts <- read.table('counts.txt', header=TRUE, sep='\t')

# 1. Normalize and compute logFC
norm_counts <- ccr.NormfoldChanges(filename='counts.txt', min_reads=30,
                                     EXPname='my_screen',
                                     libraryAnnotation=KY_Library_v1.0)

# 2. Compute genome-sorted sgRNA fold changes
gw_log_fc <- ccr.logFCs2chromPos(norm_counts$logFCs,
                                   KY_Library_v1.0)

# 3. Apply CRISPRcleanR correction
corrected <- ccr.GWclean(gw_log_fc, display=TRUE, label='my_screen')
# Output: corrected$corrected_logFCs and corrected$segments
# corrected_logFCs can replace LFCs downstream

# 4. Re-derive corrected counts for downstream MAGeCK
corrected_counts <- ccr.correctCounts('my_screen',
                                        norm_counts$norm_counts,
                                        corrected,
                                        KY_Library_v1.0,
                                        OutDir='./')
```

**Key parameter:** `min_reads=30` is the lower-count threshold for inclusion. This must match the library-coverage strategy; too high removes legitimate guides, too low keeps noisy guides.

**Output:** Pre-corrected LFCs and counts that can be fed into MAGeCK / BAGEL2 / drugZ as if they were the original screen data. The correction is independent of CN profile (unsupervised) and works on cell lines without matched WGS.


## Reconciliation: When CN Correction Fails

If post-CRISPRcleanR or post-Chronos the CN-LFC Spearman is still significantly negative, the correction is incomplete. Possible causes:

1. **Insufficient CN resolution:** A specific 4-copy region went undetected. Refine CN profile with deeper WGS.
2. **CRISPRcleanR position-based correction missed it:** The amplification is small relative to the segmentation algorithm's resolution. Use Chronos with matched CN profile.
3. **Genomic rearrangement creates a "ghost" amplification:** A complex rearrangement appears as normal CN but Cas9 cuts at multiple sites due to translocation breakpoints. Combine WGS structural variants with the analysis.
4. **Cell line has an unusually strong cut-toxicity response:** The artifact may persist; use CRISPRi screens for that line.
5. **`alternate_CN` fits one CN-effect curve across the whole panel:** it is a fitted average, not a per-line correction, so a single extreme-CN line can retain residual signal even after correction. On a real 3-line CN dose-response panel (CN=2/8/15), the CN=15 line's amplicon genes moved from -1.92 to only -1.20 -- still essential-looking -- while the essentials stayed untouched. Re-run `detect_cn_bias` per cell line after `alternate_CN`; compare those results with a pooled summary, but do not assume either view will always conceal the other.

## Apply CN Correction to Pipeline

**Workflow:**

```
1. mageck count (raw counts)
2. screen-qc verification
3. CN diagnostic: Spearman of LFC vs CN (if CN profile available)
4. If bias detected:
   a. CRISPRcleanR (pre-hoc) -> corrected counts -> MAGeCK / BAGEL2 / drugZ
   OR
   b. Chronos (joint model with CN profile) -> gene effects directly
5. Re-diagnose: Spearman of CORRECTED LFC vs CN should be near zero
6. Hit calling
```

For DepMap-style large panels:
```
Chronos handles batch + CN + screen quality in one step; no pre-correction needed.
```

For Project Score-style panel (Behan 2019):
```
CRISPRcleanR was used historically; cross-check with Chronos when CN profile available.
```

## Failure Modes

### CRISPRcleanR removes legitimate essential signal

**Trigger:** A genuine essential gene happens to lie in a region with adjacent uncorrected non-essential signal; the segment-based correction includes the essential.
**Mechanism:** CRISPRcleanR's `ccr.GWclean()` segments sgRNAs by position; segments containing multiple genes with directional consistency are corrected as a unit.
**Symptom:** A known essential drops out of post-correction hit list.
**Fix:** Inspect segments manually; if a known essential was within a corrected segment, investigate. Cross-check with non-CN-corrected MAGeCK + BAGEL2 to see if essential was a hit pre-correction.

### Chronos refuses to CN-correct fewer than 3 cell lines

**Trigger:** One or two cell lines, with or without a matched CN profile.
**Mechanism:** `alternate_CN` fits the CN-to-gene-effect curve across cell lines, so it needs a panel; its own code raises below 3.
**Symptom:** `RuntimeError: Correct for CN should not be used with fewer than 3 cell lines. Consider preprocessing with CRISPRCleanR` -- raised after a full training run, so budget for it before you spend the epochs.
**Fix:** CN-correct the counts with CRISPRcleanR pre-hoc, then hit-call; Chronos training itself works on a single line (with `negative_control_sgrnas`) and its uncorrected gene effects are still useful for QC.

### Spearman ρ still negative after CRISPRcleanR

**Trigger:** Amplification is too small or complex for the segment-based approach.
**Mechanism:** CRISPRcleanR detects systematic spatial patterns; isolated 4-copy regions can slip through.
**Symptom:** Post-correction Spearman ρ -0.05 to -0.10 between LFC and CN.
**Fix:** Refine CN profile (deeper WGS); apply Chronos with matched CN as alternative; or supplement with focal-amplification-aware methods.

### Cell line lacks matched CN profile

**Trigger:** Newly characterized line or rare patient-derived line; WGS not done.
**Mechanism:** CRISPRcleanR needs no CN profile; Chronos trains without one too, but its `alternate_CN` correction cannot run without it.
**Symptom:** Chronos gene effects still carry the artifact; CRISPRcleanR corrects unsupervised.
**Fix:** Run SNP-array (cheap, fast) or low-coverage WGS to obtain CN profile; in interim, use CRISPRcleanR unsupervised mode.

### CN amplification at non-coding region drives apparent essentiality

**Trigger:** Amplification at a gene-poor region; sgRNAs at edge genes get artifactually depleted.
**Mechanism:** Even non-essential genes adjacent to amplifications are depleted because the Cas9 cuts are at the amplified loci.
**Symptom:** Non-essential genes near amplification show LFC <0.
**Fix:** Inspect chromosomal position of "essential" hits; flag genes within 100 kb of known amplifications for orthogonal validation. This is the classic Aguirre 2016 observation.

## CRISPRi/a Alternative

**For variant-function or non-cancer-line essentiality screens, switching to CRISPRi (catalytically dead dCas9-KRAB) avoids the artifact entirely.** No DNA double-strand breaks = no DNA-damage G2 arrest = no copy-number-driven depletion.

| Approach | CN artifact | When to use |
|----------|--------------|-------------|
| Cas9 KO | YES; requires correction | Loss-of-function essentiality, traditional screens |
| CRISPRi | NO | Cancer lines with focal amps; knockdown of cuttable-toxic genes |
| CRISPRa | NO | Gain-of-function; activation screens |
| Base editing | Reduced (single-strand nick) | Variant function |
| Prime editing | Reduced | Precise edits |

See bio-crispr-screens-library-design for CRISPRi (Dolcetto) and CRISPRa (Calabrese) library options.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Spearman ρ (CN vs LFC) | <-0.10 -> bias present | Operational convention |
| Copies for detectable artifact | >6 | Operational convention; response scales with copy number (Aguirre 2016) |
| CRISPRcleanR `min_reads` | 30 (default) | Iorio 2018; lower thresholds in low-coverage screens |
| Chronos gene-effect threshold for "essential" | <-1 (cancer line) | DepMap convention |
| Chronos gene-probability for "essential" | >0.5 | DepMap convention (dependency-probability cutoff) |
| Post-correction Spearman ρ | abs(ρ) <0.05 | Acceptable correction quality |
| Cell-line CN profile resolution | ≥SNP-array level | Below this, CRISPRcleanR unsupervised |
| Cell lines needed for `chronos.alternate_CN` | ≥3 | Chronos 2.3.15 raises RuntimeError below this |
| Amplified-vs-diploid LFC gap | < -0.5 -> focal bias present | Catches single amplicons the genome-wide rho misses |
| `detect_cn_bias` low-power floor | <8 amplified genes -> distrust `bias_present: False` | Real HT-29 3-gene block: p=0.208 for a genuine -0.98 gap |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `RuntimeError: Correct for CN should not be used with fewer than 3 cell lines` | `alternate_CN` needs a panel | CN-correct with CRISPRcleanR instead; Chronos can still score the screen |
| `AssertionError: ... Is your data transposed?` | readcounts passed as guides x samples | Transpose: rows = sequence_ID, columns = sgRNA |
| `ValueError: excess_variance was passed as dict without key for '<library>': {}` at `chronos.Chronos(...)` construction (not `train()`) | `negative_control_sgrnas` not supplied | Pass a dict of non-targeting/non-essential sgRNAs per library |
| `ValueError: set of negative_control_sgrnas is empty` at `chronos.Chronos(...)` construction | `negative_control_sgrnas` was supplied, but one library's list is empty | Populate every library with real non-targeting/non-essential sgRNAs; do not use an empty placeholder list |
| CRISPRcleanR removes a known essential | Segment-based over-correction | Manually inspect segments; cross-check with non-corrected |
| Spearman ρ still -0.15 after correction | Method too coarse for the amp | Refine CN profile; use Chronos |
| ERBB2 listed as essential in SK-BR-3 | Uncorrected HER2 amplification | Always apply correction before hit calling |
| CN profile missing for newly characterized line | Profile not generated | Run SNP-array / low-coverage WGS |
| Hits restricted to non-amplified regions only | Over-correction | Reduce CRISPRcleanR aggressiveness; check known biology |
