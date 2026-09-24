---
name: bio-proteomics-quantification
category: Data Analysis
description: Quantifies protein abundance from mass spectrometry using label-free (LFQ/MaxLFQ, DIA fragment-level), isobaric (TMT/iTRAQ reporter ions, MS2 vs SPS-MS3), and metabolic (SILAC) approaches, including peptide-to-protein summarization (Tukey median polish, MaxLFQ, msqrob), sample-loading and IRS cross-plex normalization, and isotopic impurity correction. Use when turning peptide/PSM/reporter signal into a protein-by-sample abundance matrix for downstream analysis. Statistical testing of that matrix is differential-abundance; DIA quant mechanics and DIA-NN runs are dia-analysis; reading search-engine outputs is data-import; razor/shared-peptide group assignment is protein-inference.
tool_type: mixed
primary_tool: MSstats
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: MSstats 4.14.2, MSnbase 2.32.0, iq 2.0.1, numpy 2.5.3, pandas 3.0.5 (checked 2026-09-21)

Install: `pip install numpy pandas`; in R `BiocManager::install(c('MSstats', 'MSnbase'))` and `install.packages('iq')`.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Protein Quantification -- Reconstructing Protein Abundance from Ions Whose Physical Origin Dictates the Irreducible Error

**"Quantify proteins from my mass spec data"** -> Reconstruct a protein-by-sample abundance matrix from peptide/reporter ion signals, choosing a summarizer and normalizer that match where the signal physically came from -- because the measurement's physical origin sets an error that no normalization can remove.
- R: `MSstats::dataProcess()` (Tukey median polish summarization) for label-free feature-to-protein
- R: `iq::maxLFQ()` for MaxLFQ maximal peptide-ratio least-squares (normalize runs first; it does no delayed normalization)
- R: `MSnbase::quantify(reporters=TMT10)` + `purityCorrect()` for isobaric reporter extraction
- R/Python: sample-loading + IRS scaling to bridge TMT plexes; per-sample median centering for LFQ

Scope: this skill OWNS converting peptide/PSM/reporter signal into a normalized protein abundance matrix (LFQ, TMT/iTRAQ, SILAC; summarization; normalization; IRS). Statistical testing of the matrix -> differential-abundance. DIA quant mechanics and DIA-NN execution -> dia-analysis. Parsing MaxQuant/DIA-NN outputs -> data-import. Razor/shared-peptide group assignment -> protein-inference. OUT OF SCOPE: missing-value imputation and the downshift false-positive trap (modeled in differential-abundance), and absolute copy-number calibration beyond a one-line pointer. Research quantification must not be used to classify individual patients or choose treatment; route those questions to validated clinical assays.

## The Single Most Important Modern Insight

1. **Every quant method answers "where does the signal physically come from?" differently, and that physical origin dictates the error structure no normalization can remove.** LFQ measures MS1 precursor area (or DIA fragment area) in SEPARATE runs -> the irreducible error is run-to-run variation plus stochastic, left-censored (MNAR) missingness. Isobaric TMT/iTRAQ measures low-m/z reporter ions from CO-ISOLATED, co-eluting peptides in ONE spectrum -> the irreducible error is RATIO COMPRESSION toward 1:1, a PHYSICAL co-isolation effect (interloper reporters add roughly equally to every channel), attacked at the instrument by SPS-MS3 and never fully undone in software. SILAC measures a heavy/light MS1 pair in the SAME scan -> lowest per-ratio variance, but its irreducible vulnerabilities are incomplete labeling and Arg->Pro label scrambling, which bias every ratio and cannot be corrected post hoc because channels are combined before any MS (a correctly mixed sample carries no inherent mixing error). Teach the signal origin and every threshold and failure mode below follows from it.

2. **The peptide-to-protein SUMMARIZATION choice is the highest-leverage decision in the pipeline, and it is invisible in the output.** In log space a peptide intensity is protein abundance + a peptide effect (ionization efficiency, flyability, missed cleavages, modifications) that spans orders of magnitude and is partly context-dependent, plus a run effect. Sum is dominated by the highest-flying peptide (dropout collapses it -> a fold change driven by detectability, not biology); mean is unbiased only if the detected peptide SET is identical across runs (it is not under MNAR); median discards relative-intensity information. Benchmarks confirm the quantification method is a dominant driver of which proteins are called differential (Lin 2022). Report the summarizer as prominently as the test, and run a sensitivity analysis across >=2 summarizers -- that is where the answer is most likely to move.

## Tool Taxonomy

| Tool / method | Citation | Mechanism / role | When |
|---|---|---|---|
| MaxLFQ | Cox 2014 | delayed normalization + maximal shared-peptide log-ratio least-squares; peptide scale cancels in pairwise ratios | label-free DDA/DIA relative quant across many samples |
| `iq::maxLFQ()` | Cox 2014; Pham 2020 | R implementation of MaxLFQ ratio extraction (no delayed normalization; give it run-normalized log2 input); call it, do NOT reimplement | running MaxLFQ outside MaxQuant/DIA-NN |
| MSstats Tukey median polish | Choi 2014 | iteratively subtract peptide-medians + run-medians in log space; column effects = per-sample abundance; 50% breakdown | robust label-free default summarizer |
| msqrob (peptide-level) | Sticker 2020; Goeminne 2016 | treats the peptide effect as a covariate not noise; ridge + empirical Bayes + Huber | accuracy-critical small-n, unbalanced coverage (route OUT to differential-abundance) |
| iBAQ | -- | sum(peptide intensities) / number of theoretically observable tryptic peptides | rank / order-of-magnitude within-sample abundance |
| Top3 / Hi3 | Silva 2006 | sum/avg of top-3 peptide intensities, calibrate with one spiked standard | absolute amount, ~2 orders linear |
| Proteomic ruler | Wisniewski 2014 | histone signal as internal molar reference, no spike-in | absolute copies/cell without standards |
| Spectral counting / NSAF | Zybailov 2006 | count PSMs per protein, divide by protein LENGTH then total | largely OBSOLETE; niche AP-MS only |
| TMT/iTRAQ reporter | Ting 2011; McAlister 2014 | isobaric tag; reporter ratios at MS2 or SPS-MS3 | high multiplexing, zero run-to-run variation within a plex |
| SILAC | Ong 2002 | heavy/light precursor pair co-elute in the SAME MS1 scan | lowest-variance ratios, cell culture that can be labeled |
| DIA fragment-level | Demichev 2020 | MaxLFQ at the FRAGMENT level then roll up (route OUT to dia-analysis) | low-missingness label-free cohorts |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|---|---|---|
| Label-free DDA, MaxQuant evidence.txt | `MSstats::dataProcess` (TMP) | robust summarization, the workhorse; censored (AFT) imputation only with `MBimpute = TRUE` |
| Label-free, need MaxLFQ outside MaxQuant | `iq::maxLFQ()` | the real algorithm; median centering is NOT MaxLFQ |
| Label-free DIA matrix | -> dia-analysis | DIA-NN MaxLFQ at fragment level is owned there |
| TMT, accuracy critical | SPS-MS3 acquisition + reporter extraction | co-isolation rejected at the instrument; software cannot fully undo compression |
| TMT, single plex only | MS2 reporters + sample-loading normalization (`references/tmt_isobaric.md`) | within-plex ratios are stable |
| TMT, multiple plexes | sample-loading THEN IRS bridge (Plubell 2017; `references/tmt_isobaric.md`) | absolute reporter intensities are NOT comparable across runs without a reference channel |
| SILAC ratios | verify labeling efficiency + Arg->Pro first (`references/silac.md`) | unchecked, both bias every ratio invisibly |
| Absolute copy number | proteomic ruler or Top3 + standard | iBAQ is within-sample rank only |
| AP-MS / affinity-enrichment pulldown | do NOT median/SL/IRS-normalize; control subtraction (SAINT/CompPASS/CRAPome; `references/affinity_enrichment.md`) | an enrichment is not a balanced proteome; data-internal normalization erases the bait signal |
| Which summarizer? | run >=2 (TMP and MaxLFQ) and compare | this is the highest-leverage, invisible choice |

Default when uncertain: label-free DDA -> `MSstats::dataProcess` with `summaryMethod='TMP'`, `normalization='equalizeMedians'`; report the summarizer alongside results and sanity-check against `iq::maxLFQ()`.

## Label-Free Summarization and Normalization

### Summarize peptides to proteins with MSstats

**Goal:** Turn MaxQuant feature-level evidence into a normalized protein-level abundance matrix.

**Approach:** Read the MaxQuant tables with quoting off (protein names contain apostrophes), reformat to MSstats input, then `dataProcess` applies median equalization and Tukey median polish (robust to outlier peptides, 50% breakdown). With `MBimpute = FALSE` there is no censored-value model; `MBimpute = TRUE` is MSstats' AFT imputation of censored features (the only censoring route in MSstats; `groupComparison` has none). Otherwise model missingness downstream (proDA/msqrob2 in differential-abundance).

```bash
Rscript scripts/msstats_summarize.R evidence.txt proteinGroups.txt annotation.csv protein_abundance.csv   # optional 5th arg: TRUE = MBimpute
```
The script reads the tables with `quote = ''` and checks the row count, then runs `MaxQtoMSstatsFormat` and `dataProcess(normalization = 'equalizeMedians', summaryMethod = 'TMP', censoredInt = 'NA', MBimpute = FALSE)`; on/off proteins come out of `groupComparison` as log2FC `-Inf` / issue `oneConditionMissing` when `MBimpute = FALSE`. It writes `ProteinLevelData`.

### Run the real MaxLFQ (not median centering)

**Goal:** Produce MaxLFQ protein intensities from a peptide quant matrix.

**Approach:** Call `iq::maxLFQ()`, which implements the Cox 2014 maximal peptide-ratio least-squares step. It does NOT perform delayed normalization, so run-level loading offsets pass straight into the estimates: median-normalize each run's peptide log2 intensities first. `maxLFQ()` takes ONE protein's matrix, so a whole table goes through `preprocess` -> `create_protein_list` -> `create_protein_table`. Per-sample median centering of the protein matrix shares only the name and silently gives a different answer.

```bash
Rscript scripts/maxlfq_iq.R peptide_long.csv protein_maxlfq.csv   # columns protein, ion, run, intensity (RAW)
```
Runs `preprocess(median_normalization = TRUE)` -> `create_protein_list` -> `create_protein_table(method = 'maxLFQ')` and writes the proteins x samples `$estimate`. `maxLFQ()` solves per-protein least squares only within CONNECTED sample sets; a non-empty `$annotation` marks proteins whose samples are NOT on one common scale, and those are listed in `<out>.disconnected.txt`.

One protein at a time, in R:
```r
library(iq)
# protein_list from create_protein_list(); rows = peptide ions, columns = samples, values = RUN-NORMALIZED log2 intensities
result <- maxLFQ(protein_list[[1]])
# $estimate is an UNNAMED vector in the input column order; name it or samples silently transpose
protein_estimate <- setNames(result$estimate, colnames(protein_list[[1]]))
```

### Median-center label-free intensities (a normalizer, not a summarizer)

**Goal:** Correct per-sample loading differences before testing.

**Approach:** Subtract each sample's median log2 intensity (corrects LOCATION only and cannot manufacture variance). It assumes most proteins are unchanged and changes are roughly symmetric: an up-heavy quantified set shifts every fold change, so check the median fold change of presumed-null proteins or normalize on a stable subset. Median centering normalizes; it does NOT summarize peptides to proteins.

```python
import numpy as np
import pandas as pd

log_int = np.log2(intensities.replace(0, np.nan))    # MaxQuant writes 0 for missing; log2(0) = -inf
sample_medians = log_int.median(axis=0)
normalized = log_int - sample_medians + sample_medians.median()
```

Runnable with asserted output: `examples/lfq_normalization.py` (median centering, sample-loading + IRS, the SILAC pilot checks and the AP-MS scorer on seeded data).

## Reference Files

Read only the file the request needs; the LFQ/MSstats/MaxLFQ core, decision tree, failure modes, thresholds and Common Errors stay here.

| File | Read when |
|---|---|
| `references/tmt_isobaric.md` | TMT/iTRAQ reporter extraction, lot-specific impurity correction (incl. TMTpro CoA route), sample-loading + IRS bridging of multiple plexes |
| `references/silac.md` | SILAC: labeling-efficiency and Arg->Pro checks on a heavy-only pilot, log2 H/L ratios that keep on/off proteins |
| `references/affinity_enrichment.md` | AP-MS / affinity-enrichment / proximity-labeling pulldown: scoring prey against negative-control IPs |

## Per-Method Failure Modes

### MaxLFQ reimplemented as median centering
**Trigger:** A homebrew function named `maxlfq` that only subtracts per-sample medians.
**Mechanism:** Real MaxLFQ is delayed normalization plus a maximal peptide-ratio least-squares solve; median centering shares only the name.
**Symptom:** Plausible-looking but systematically different intensities; unbalanced peptide sets handled wrongly.
**Fix:** Call `iq::maxLFQ()` on run-normalized input, DIA-NN, or MaxQuant.

### TMT ratio compression
**Trigger:** MS2-only reporter quant on a complex sample.
**Mechanism:** Co-isolated interloper peptides add reporters roughly equally to every channel, pulling large true ratios toward 1:1; PHYSICAL, not removable by normalization.
**Symptom:** "Nothing is significant"; attenuated fold changes, inflated false negatives.
**Fix:** SPS-MS3 acquisition, narrower isolation windows, FAIMS/ion mobility, or complement-reporter methods; a PIF filter helps but MS1 purity underestimates true interference (Savitski 2013).

### TMT plexes concatenated without IRS
**Trigger:** Stacking reporter intensities from multiple plexes directly.
**Mechanism:** Absolute reporter intensities for one protein differ 2-5x across runs from random elution-profile sampling, unrelated to abundance.
**Symptom:** Plex appears as the dominant axis of variation; spurious cross-plex differences.
**Fix:** Include a pooled reference channel in EVERY plex; apply sample-loading then IRS (Plubell 2017).

### Isotopic impurity matrix mis-applied
**Trigger:** Using the default/example impurity matrix, the wrong lot, or a transposed/mis-ordered (127N vs 127C) matrix.
**Mechanism:** A few percent of each channel bleeds to +/-1 Da neighbors; wrong values mis-subtract, negative corrected intensities get clipped. In TMT10 order (126, 127N, 127C, ...) N and C channels interleave, so the 13C +1 Da neighbour of 126 is 127C, two positions away -- a CoA sheet laid out as adjacent columns is misplaced.
**Symptom:** Adjacent channels silently biased; extreme contrasts placed in adjacent channels confounded.
**Fix:** Use the lot-specific Certificate of Analysis values; randomize channel-to-condition assignment.

### SILAC Arg->Pro conversion / incomplete labeling
**Trigger:** Pro-containing peptides or a labeling efficiency below ~95%.
**Mechanism:** Cells convert heavy Arg to heavy Pro (+6 Da), splitting Pro-peptide signal and underestimating the heavy channel; residual light masquerades as down-regulation.
**Symptom:** Ratios biased toward light, worse for Pro-rich proteins; invisible without a check.
**Fix:** Proline supplementation, measure conversion per cell line, verify >=95% incorporation on a heavy-only pilot -- `silac_labeling_efficiency` and `arg_to_pro_shift` in `references/silac.md` compute both.

### SILAC on/off proteins discarded as NaN
**Trigger:** Returning NaN whenever either channel is zero.
**Mechanism:** A protein present only in heavy (or only light) is the interesting biology, thrown away.
**Symptom:** Largest true changes silently dropped before analysis.
**Fix:** Flag present-in-one-channel cases (H-only / L-only) in a separate column and report them as a presence table; do not put +/-Inf into the matrix handed to limma. Route honest absence handling to differential-abundance.

### Spectral counting / NSAF used for fold changes
**Trigger:** Reaching for PSM counts for quantitative comparison.
**Mechanism:** Count statistics are catastrophic at low abundance and saturate; dynamic exclusion deliberately breaks count-abundance proportionality. NSAF divides intensity by protein LENGTH then total -- a count divided by total spectra is not NSAF.
**Symptom:** Noisy, biased estimates; the wrong normalization labeled NSAF.
**Fix:** Use MS1/MS2 intensity (LFQ/DIA); keep spectral counting as historical context only.

### AP-MS / enrichment pulldown normalized as a balanced proteome
**Trigger:** Median/sample-loading/IRS normalization applied to an affinity-purification or biotin-enrichment pulldown.
**Mechanism:** Data-internal normalization assumes most signal is an unchanging background; a successful pulldown is deliberately non-representative (bait plus a few interactors over background), so equalizing medians or loading rescales away the enrichment being measured.
**Symptom:** Real interactors flattened toward background; bait abundance dominates the axis of variation.
**Fix:** Do not data-internal-normalize an enrichment; score against negative-control pulldowns with `score_vs_control_ips` in `references/affinity_enrichment.md` (SAINTexpress/CompPASS/CRAPome for a probability) or normalize to bait abundance.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| MaxLFQ min. ratio count = 2 | Cox 2014 (default) | a single shared peptide gives a ratio with no outlier rejection; >=2 lets the median start rejecting interference; setting 1 admits unguarded single-peptide ratios |
| PIF >= 0.75 | community filter | rejects spectra with too much interloper signal; but MS1 purity UNDERESTIMATES true reporter interference (Savitski 2013), so PIF ~0.9 can still be compressed |
| TMT N/C reporter spacing = 6.3 mDa | Thompson 2019 | 13C-vs-15N mass defect; needs high-res MS2 (>=30-50k) to resolve N from C channels |
| Reporter match tolerance ~0.002-0.003 Da | -- | tight enough to separate 6.3 mDa N/C channels at high resolution |
| SILAC labeling efficiency >= 95% (target 97-98%) | -- | residual light contaminates the heavy channel -> false down-regulation; needs ~5-6 doublings |
| Min peptides per protein for quant >= 2 | -- | single-peptide (one-hit) proteins are quant-unreliable |
| Top3 uses exactly the top 3 peptides | Silva 2006 | most intense peptides are most reproducibly detected, closest to uniform per-mole response |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `readMSnSet` does not extract reporters | it reads an already-quantified text matrix | use `readMSData(mode='onDisk')` then `quantify(reporters=TMT10, method='max')` |
| `log2(0) = -inf` in the matrix | MaxQuant writes 0 for "not quantified" | replace 0 -> NaN before any transform |
| Reading `Intensity` when ratios needed | `Intensity` is raw, not normalized; `iBAQ` is within-sample only | use `LFQ intensity` for between-sample LFQ ratios (see data-import) |
| Median centering called MaxLFQ | homebrew shares only the name | call `iq::maxLFQ()` / DIA-NN / MaxQuant |
| `maxLFQ()` returns 8 numbers with no sample names; columns end up scrambled | `$estimate` is an unnamed vector in the input column order, and `maxLFQ()` handles ONE protein | `setNames(result$estimate, colnames(m))`; for a whole table use `preprocess` -> `create_protein_list` -> `create_protein_table` |
| On/off proteins get `log2FC = -Inf`, `issue = oneConditionMissing`; no censored handling | `MBimpute=FALSE` means no censoring model; MSstats' AFT exists only in `dataProcess(MBimpute=TRUE)` | keep `MBimpute=TRUE` for the MSstats AFT route, or keep `FALSE` and model missingness with proDA/msqrob2 (differential-abundance) |
| `EOF within quoted string` and far fewer rows than the file | default `read.table` quoting on MaxQuant tables with apostrophes | `read.table(..., quote = '', comment.char = '')` or `data.table::fread`; check the row count |
| Script hangs at `makeImpuritiesMatrix` under Rscript | default `edit = TRUE` opens an editor | `makeImpuritiesMatrix(x = 10, edit = FALSE)` or `filename = 'lot_coa.csv', edit = FALSE` |
| `makeImpuritiesMatrix(x = 16)`: `length of 'dimnames' [1] not equal to array extent` | MSnbase templates exist only for x = 4, 6, 8, 10; TMTpro falls through to an unnamed `diag(x)` | supply the lot CoA: `makeImpuritiesMatrix(filename = 'lot_coa.csv', edit = FALSE)` with 16 offset columns; quantify with `reporters = TMT16` (there is no `TMT18`) |
| Whole plex rows become Inf/NaN after IRS | reference channel 0 or missing for that protein | mask references `<= 0`/NaN and report the unbridged proteins |
| `eBayes`: `prior.weights contain NA values` on SILAC ratios | +/-Inf on/off codes inside the ratio matrix | NaN ratio plus a presence flag column; drop rows with no finite ratio before `lmFit` |
| Cross-plex TMT comparison is invalid | no reference channel / no IRS | add a pooled reference channel per plex, apply SL then IRS |
| Every unchanged protein sits below 0 in log2 H/L | incomplete heavy incorporation: the unlabeled (1 - e) fraction is counted as light | measure e on a heavy-only pilot (`silac_labeling_efficiency`); a 1:1 mix then reads log2 H/L = log2(e / (2 - e)), -0.20 at e = 0.93 |
| Pro-rich proteins look systematically down in SILAC | Arg->Pro conversion, once per proline | `arg_to_pro_shift`: a slope in log2 H/L per proline (a flat offset is labeling, not scrambling); supplement proline |
| AP-MS "top N by fold change over the input lysate" is all ribosome/chaperone/tubulin | the input is not a negative control; bead background is enriched in the IP and diluted in the input | score against control IPs (`score_vs_control_ips`) and require detection in every bait replicate |
| MSnbase deprecation warnings | MSnbase is in maintenance mode | current pipelines use Spectra + QFeatures (`readQFeatures`, `aggregateFeatures`) |

## References

- Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. 2014. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. *Mol Cell Proteomics* 13(9):2513-2526.
- Silva JC, Gorenstein MV, Li GZ, Vissers JPC, Geromanos SJ. 2006. Absolute quantification of proteins by LCMSE: a virtue of parallel MS acquisition. *Mol Cell Proteomics* 5(1):144-156.
- Wisniewski JR, Hein MY, Cox J, Mann M. 2014. A "proteomic ruler" for protein copy number and concentration estimation without spike-in standards. *Mol Cell Proteomics* 13(12):3497-3506.
- Zybailov B, Mosley AL, Sardiu ME, et al. 2006. Statistical analysis of membrane proteome expression changes in Saccharomyces cerevisiae. *J Proteome Res* 5(9):2339-2347.
- Ong SE, Blagoev B, Kratchmarova I, et al. 2002. Stable isotope labeling by amino acids in cell culture, SILAC, as a simple and accurate approach to expression proteomics. *Mol Cell Proteomics* 1(5):376-386.
- Ting L, Rad R, Gygi SP, Haas W. 2011. MS3 eliminates ratio distortion in isobaric multiplexed quantitative proteomics. *Nat Methods* 8(11):937-940.
- McAlister GC, Nusinow DP, Jedrychowski MP, Wuhr M, et al. 2014. MultiNotch MS3 enables accurate, sensitive, and multiplexed detection of differential expression across cancer cell line proteomes. *Anal Chem* 86(14):7150-7158.
- Savitski MM, Mathieson T, Zinn N, et al. 2013. Measuring and managing ratio compression for accurate iTRAQ/TMT quantification. *J Proteome Res* 12(8):3586-3598.
- Thompson A, Wolmer N, Koncarevic S, et al. 2019. TMTpro: design, synthesis, and initial evaluation of a proline-based isobaric 16-plex tandem mass tag reagent set. *Anal Chem* 91(24):15941-15950.
- Plubell DL, Wilmarth PA, Zhao Y, et al. 2017. Extended multiplexing of tandem mass tags (TMT) labeling reveals age and high fat diet specific proteome changes in mouse epididymal adipose tissue. *Mol Cell Proteomics* 16(5):873-890.
- Choi M, Chang CY, Clough T, et al. 2014. MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments. *Bioinformatics* 30(17):2524-2526.
- Goeminne LJE, Gevaert K, Clement L. 2016. Peptide-level robust ridge regression improves estimation, sensitivity, and specificity in data-dependent quantitative label-free shotgun proteomics. *Mol Cell Proteomics* 15(2):657-668.
- Sticker A, Goeminne L, Martens L, Clement L. 2020. Robust summarization and inference in proteome-wide label-free quantification. *Mol Cell Proteomics* 19(7):1209-1219.
- Demichev V, Messner CB, Vernardis SI, Lilley KS, Ralser M. 2020. DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. *Nat Methods* 17(1):41-44.
- Pham TV, Henneman AA, Jimenez CR. 2020. iq: an R package to estimate relative protein abundances from ion quantification in DIA-MS-based proteomics. *Bioinformatics* 36(8):2611-2613.
- Ong SE, Mann M. 2006. A practical recipe for stable isotope labeling by amino acids in cell culture (SILAC). *Nat Protoc* 1(6):2650-2660. (Incorporation check and Arg->Pro.)
- Sowa ME, Bennett EJ, Gygi SP, Harper JW. 2009. Defining the human deubiquitinating enzyme interaction landscape. *Cell* 138(2):389-403. (CompPASS.)
- Teo G, Liu G, Zhang J, Nesvizhskii AI, Gingras AC, Choi H. 2014. SAINTexpress: improvements and additional features in Significance Analysis of INTeractome software. *J Proteomics* 100:37-43.
- Mellacheruvu D, Wright Z, Couzens AL, et al. 2013. The CRAPome: a contaminant repository for affinity purification mass spectrometry data. *Nat Methods* 10(8):730-736.
- Lin MH, Wu PS, Wong TH, Lin IY, Lin J, Cox J, Yu SH. 2022. Benchmarking differential expression, imputation and quantification methods for proteomics data. *Brief Bioinform* 23(3):bbac138.

## Related Skills

- data-import - Parse MaxQuant/DIA-NN outputs and pick the right intensity column before quantifying
- protein-inference - Razor/shared-peptide group assignment that determines which protein a peptide counts toward
- differential-abundance - Statistical testing, missing-value modeling, and the downshift false-positive trap
- proteomics-qc - CV, correlation, and PCA checks that confirm normalization worked
- dia-analysis - DIA fragment-level MaxLFQ and DIA-NN execution
- differential-expression/de-results - Shared empirical-Bayes and FDR conventions for expression matrices
- data-visualization/heatmaps-clustering - Visualize the normalized abundance matrix
- workflows/proteomics-pipeline - End-to-end pipeline that calls this skill for the quant step
