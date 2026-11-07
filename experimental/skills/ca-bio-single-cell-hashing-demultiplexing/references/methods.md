# Methods and limitations

Source tool versions and synthetic observations are unverified here. Check completeness and errors together; a partial output file alone does not prove a successful run. Multiplet proportions also depend on pooling balance.

## Governing principle

In cell or nucleus hashing, each sample is labeled before pooling with a unique oligo-tagged reagent: a barcoded antibody (CITE-seq HTO), a lipid- or cholesterol-modified oligo (MULTI-seq), or a CellPlex CMO. A true singlet's HTO count vector is dominated by ONE tag standing well above a background of ambient and spillover counts, so sample assignment reduces to one question per cell: which tag, if any, exceeds that cell's background. Cross-sample doublets fall out directly - two tags both high - which is the decisive advantage of hashing over expression-only doublet detection.

Hashtag demux, genetic demux, and expression-doublet detection answer three DIFFERENT questions and should be combined, not substituted. Hashtag and genetic demux both assign samples and catch cross-sample doublets, but expression-doublet detection (single-cell/doublet-detection) catches within-sample and homotypic doublets that hashing and genetics are blind to, because two cells from the same sample carry the same tag and the same genotype. Conversely, expression methods miss cross-sample doublets when the two samples are transcriptionally similar. The cross-sample doublet rate also calibrates the expected TOTAL doublet rate: with two pooled samples within- and cross-sample doublets are equally frequent, while with k samples cross-sample doublets dominate and within-sample ones fall to about 1/k of all doublets, so hashing still misses that within-sample fraction (expression-doublet detection stays necessary) and a hashing doublet rate far below the expression-doublet rate is a red flag that staining or thresholds are off.

The hard part is the background, not the dominant tag. Ambient HTO from lysed cells, spillover between tags, staining failure, and batch differences in tag-capture efficiency all inflate the background and grow the "Negative" pile (real cells whose true tag never cleared background, distinct from true empty droplets removed upstream). Nucleus hashing is harder than whole-cell because tag capture is lower. Methods differ mainly in how they model that background.

## Choosing a demultiplexing modality

| Modality | Needs | Cross-sample doublets | Cannot do | Tools |
|----------|-------|-----------------------|-----------|-------|
| Hashtag/HTO (this skill) | HTO/lipid/CMO library at pooling | Yes (two tags high) | Nothing without a hashing library; sensitive to staining and ambient | HTODemux, MULTIseqDemux, hashsolo, demuxEM, GMM-Demux, demuxmix |
| Genetic (natural SNPs) | >=2 distinct genotypes, no hashing | Yes | Cannot separate same-donor samples (identical genotype) | demuxlet, freemuxlet, souporcell, vireo |
| Expression doublet (orthogonal) | Just the GEX matrix | No (misses when samples similar) | Catches within-sample/homotypic doublets the other two miss | scDblFinder, Scrublet (single-cell/doublet-detection) |

Genetic demux is the fallback when no hashing was done but samples come from different donors; it cannot resolve multiple samples from one donor, which hashing can. Pair whichever sample-assignment method applies with expression-doublet detection for the doublets it cannot see.

## Choosing a hashtag caller

| Method | Model | Use when | Fails when |
|--------|-------|----------|------------|
| HTODemux (Seurat) | k-medoids cluster per HTO + negative-distribution quantile | Standard antibody HTO, clean bimodal staining, Seurat workflow | Weak/low-depth staining or heavy ambient; clustering unstable on near-zero HTOs |
| MULTIseqDemux (Seurat) | Per-HTO KDE, threshold between maxima, quantile sweep | MULTI-seq lipid/cholesterol tags; want autoThresh to optimize the quantile | Few cells; unimodal density when one tag dominates |
| hashsolo (scanpy/solo) | Bayesian over negative/singlet/doublet | Few hashtags (works at 2), many negatives, scanpy-native pipeline | Very low signal; priors mis-set for the actual doublet rate |
| demuxEM (pegasus) | EM with background estimated from empty droplets | High ambient; nucleus hashing; raw matrix with empties available | Empty droplets filtered out before calling; very sparse signal |
| GMM-Demux (CLI) | Gaussian mixture on normalized HTO, explicit MSM multiplets | Want explicit multiplet accounting or experiment planning | Non-Gaussian background; poor per-tag separation |
| demuxmix (R) | Negative-binomial regression mixture, optional RNA covariate | Bad/variable staining, batch tag-efficiency differences | Very few cells per mixture component |

The EM and regression methods (demuxEM, demuxmix) model the background explicitly and are the robust choice when staining is marginal, and they handle a two-tag pool as readily as a many-tag one (their failure mode is too few cells per mixture component, not too few tags), so weak staining even at two tags routes to demuxmix or demuxEM rather than hashsolo; HTODemux and MULTIseqDemux are fast defaults for clean data; hashsolo handles few hashes and many negatives. When callers disagree, run a consensus (cellhashR wraps several callers) and verify current best practice against installed docs before trusting any single call.

## Model ambient background explicitly (Python / R)

**Goal:** Recover correct calls when ambient HTO or weak staining inflates the background.

**Approach:** demuxEM estimates the background from empty droplets before assigning signal; demuxmix fits a negative-binomial regression mixture using the number of detected genes as a covariate, both of which are more robust than a fixed quantile. demuxEM's `pegasusio` dependency is a compiled C extension that does not build on Windows (its `getline` call has no MSVC equivalent) - verify the pattern below against `pegasus`'s installed docstrings before trusting it on a machine where it hasn't been introspected; on Linux/Mac, `import pegasus; help(pegasus.demultiplex)` first.

```python
import pegasus as pg

pg.estimate_background_probs(hashing_data)
pg.demultiplex(rna_data, hashing_data, min_signal=10.0)
rna_data.obs['demux_type'].value_counts()    # singlet / doublet / unknown
rna_data.obs['assignment']                    # sample name per cell
```

```r
library(demuxmix)

dmm <- demuxmix(as.matrix(hto_counts), rna = num_detected_genes)
calls <- dmmClassify(dmm)                      # HTO assignment + Type (singlet/multiplet/negative/uncertain)
```

`min_signal=10.0` marks cells with too little signal as unknown; lower it to rescue low-capture nucleus hashing, raise it for cleaner singlets. demuxmix's RNA covariate is what makes it robust to per-tag staining differences.

**demuxmix can hard-crash instead of degrading gracefully.** On HTO counts that are underdispersed relative to what its per-tag negative-binomial fit expects, `demuxmix()` throws an uncaught R error from `glm.nb` (`missing value where TRUE/FALSE needed`) rather than returning a result - confirmed by reproducing the crash on synthetic near-Poisson background counts. demuxmix's own warning fires first and names the cause: `Underdispersion observed for HTO "..."; consider running demuxmix with a manual initial droplet assignment using the clusterInit parameter`. On the same data that crashed, `model = 'naive'` (drop the RNA-regression covariate) completed instead of crashing - confirmed by re-running the identical input - though `dmmClassify()` may still warn `Not all models converged. Do not use the classification results.`, which must be checked before trusting the output:

```r
dmm <- tryCatch(
  demuxmix(as.matrix(hto_counts), rna = num_detected_genes),
  error = function(e) {
    message('demuxmix regression fit failed (', conditionMessage(e), '); retrying with model="naive"')
    demuxmix(as.matrix(hto_counts), model = 'naive')
  }
)
calls <- dmmClassify(dmm)   # check for a "did not converge" warning before trusting calls
```

## Demultiplex with explicit multiplet accounting: GMM-Demux (CLI)

**Goal:** Gaussian-mixture classification per tag with an explicit multi-sample-multiplet (MSM) breakdown, independent of Seurat/scanpy.

**Approach:** GMM-Demux reads a cells-by-HTO CSV or a CellRanger mtx folder, fits a 2-component Gaussian mixture per tag, and writes a per-cell `Cluster_id` (0 = negative, 1..n = each HTO singlet, remaining = every multiplet combination, decoded by the accompanying `.config` file). Checked on GMM-Demux 0.2.2.3: verified end to end on a synthetic 1600-cell, 4-tag dataset (1600/1600 global and 1360/1360 singlet sample-ID agreement against ground truth).

```bash
GMM-demux -c hto_counts.csv HTO_A,HTO_B,HTO_C,HTO_D -f gmm_out
```

`-c` reads CSV (cells x HTO columns, first column the cell barcode) instead of the default mtx format; `-f gmm_out` writes `GMM_full.csv` (`Cluster_id`, `Confidence` per cell) and `GMM_full.config` (the `Cluster_id` -> tag-combination legend) into `gmm_out/`. Two caveats found by running it: input HTO columns must be a float dtype in the CSV, not integer, or GMM-Demux's internal CLR step throws `LossySetitemError` on recent pandas - cast with `df[hto_cols] = df[hto_cols].astype(float)` before writing the CSV if you hit that. Separately, after classification GMM-Demux writes a same-sample-droplet mtx folder (`SSD_mtx/` in the working directory unless `-o <dir>` is given), and on some inputs or scipy versions that step has thrown `ValueError: unsupported data types in input` even though the `-f` report was already written correctly and the run's actual result was valid. It is not deterministic: on GMM-Demux 0.2.2.3 with scipy 1.13.1 it exited 0 on 3 of 3 synthetic datasets (1600-cell 4-tag, unequal-pooling 4-tag, weak-staining 3-tag). **Check that `GMM_full.csv` exists and is non-empty rather than trusting the exit code.**

## Threshold and parameter reference

| Parameter | Default | Rationale and when to change |
|-----------|---------|------------------------------|
| HTODemux positive.quantile | 0.99 | Quantile of the negative distribution defining "positive"; raise for stricter calls (more Negatives), lower to rescue weak staining |
| NormalizeData CLR margin | 1 (Seurat default); 2 common for HTO | margin=2 normalizes each tag across cells, correcting per-tag capture bias; pick per the staining and verify |
| MULTIseqDemux quantile / autoThresh | 0.7 / FALSE | autoThresh sweeps the quantile to maximize singlets; use when a fixed threshold over- or under-calls |
| hashsolo priors | (0.01, 0.8, 0.19) | [negative, singlet, doublet]; the doublet prior should track expected loading doublets (~0.8% per 1000 cells on 10x) |
| hashsolo number_of_noise_barcodes | `len(cell_hashing_columns) - 2` | **0 at 2 tags, 1 at 3 tags - set explicitly (e.g. 1) whenever you have 2-3 hashtags**, or the noise fit degenerates and every cell is silently called Negative |
| demuxEM min_signal | 10.0 | Cells below this signal are unknown; lower for low-capture nuclei, raise for cleaner singlets |
| demuxmix clusterInit / model | `list()` / "auto" | On an underdispersion warning or a glm.nb crash, retry with `model='naive'` or a manual `clusterInit`, per demuxmix's own warning |
| GMM-Demux -t threshold | 0.8 | Confidence threshold for a droplet to be assigned rather than left ambiguous; lower to rescue more droplets at the cost of confidence |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Huge Negative pile, few singlets | Weak staining or high ambient inflating background; quantile too strict | Lower positive.quantile / min_signal; switch to demuxEM (empty-droplet background) or demuxmix (RNA covariate) |
| Cross-sample doublet rate near zero but expression doublets high | Hashing thresholds too loose, or doublet prior too low | Tighten the quantile; raise hashsolo doublet prior; reconcile against the expected loading doublet rate |
| HTODemux errors on a zero-count cluster | Cells with all-zero HTO counts cluster together | Filter cells with no HTO counts before HTODemux; check the HTO matrix barcodes match the GEX cells |
| Calls flip with normalization choice | CLR margin=1 vs margin=2 shifts per-tag thresholds | Choose margin deliberately (margin=2 corrects tag-efficiency bias); compare both and inspect ridge plots |
| Genetic demux cannot split two samples | Both samples are the same donor (identical genotype) | Use hashtag demux; genetic methods cannot separate same-donor samples |
| Cross-sample doublets present but homotypic doublets remain | Hashing is blind to within-sample doublets | Run expression-doublet detection (single-cell/doublet-detection) in addition |
| Nucleus hashing yields mostly Negatives | Lower tag capture in nuclei than whole cells | Use demuxEM (designed for nuclei); lower min_signal; expect a larger Negative fraction |
| Looks clean (low Negatives, low doublets) but nearly all cells are one sample | Staining failure where one tag dominates all cells (mispipetted/over-concentrated antibody, or all samples got the same tag) | Sanity-check the per-tag singlet distribution against the expected pooling; one tag capturing nearly all cells means staining failed even though Negatives look low |
| One sample silently lost or contaminating while others demultiplex fine | A single antibody failed to stain, so its cells fall into Negative or misassign to the nearest-ambient tag | Check each tag has a non-trivial positive population; one near-zero tag means a failed antibody dropped or misassigned that sample |
| The rare sample in unequal pooling is under-recovered | Very unequal pooling (e.g. 80/10/10) leaves the minority tag too few positives to form a clean cluster or negative distribution | Inspect per-tag ridge plots; consider demuxmix for the minority tag, whose regression mixture is more stable on small components |
| Many cells flagged generic "Doublet" | Cells positive for 3+ tags collapsed to one label, hiding over-loading or heavy ambient | Inspect the multiplet tag-count distribution (GMM-Demux MSM); 3+ tags high is a run-quality diagnostic, not an ordinary 2-cell doublet |
| hashsolo classifies (almost) every cell "Negative" with 2-3 hashtags, no error raised | `number_of_noise_barcodes` defaults to `len(cell_hashing_columns) - 2` = 0 at 2 tags, degenerating the noise fit; only an unrelated numpy "Mean of empty slice" warning appears | Pass `number_of_noise_barcodes=1` explicitly at 2-3 tags; check the Negative fraction is not near 100% before trusting the result |
| demuxmix throws an uncaught `glm.nb` error (`missing value where TRUE/FALSE needed`) | HTO background counts are underdispersed for the per-tag negative-binomial fit; demuxmix's own warning names it before crashing | Retry with `model='naive'` (verified to complete on the same data that crashed) or a manual `clusterInit`; check `dmmClassify()`'s convergence warning either way |
| GMM-Demux exits non-zero but `GMM_full.csv` looks populated | Its SSD-mtx writer step runs after classification and can fail on some inputs/scipy versions, after the classification report was already written correctly | Inspect the error and output completeness together; a partial file does not establish success |

## References

- Stoeckius et al. 2018, Genome Biol 19:224 - Cell Hashing; barcoded antibodies for multiplexing and doublet detection.
- McGinnis et al. 2019, Nat Methods 16:619-626 - MULTI-seq; lipid- and cholesterol-tagged-oligo sample multiplexing.
- Kang et al. 2018, Nat Biotechnol 36:89-94 - demuxlet; genetic demultiplexing from natural variation.
- Heaton et al. 2020, Nat Methods 17:615-620 - souporcell; genotype clustering without reference genotypes.
- Huang et al. 2019, Genome Biol 20:273 - vireo; Bayesian genetic demultiplexing without a genotype reference.
- Bernstein et al. 2020, Cell Syst 11(1):95-101 - Solo and hashsolo; Bayesian hashing demultiplexing.
- Gaublomme et al. 2019, Nat Commun 10:2907 - demuxEM; nuclei multiplexing with background estimated from empty droplets.
- Xin et al. 2020, Genome Biol 21:188 - GMM-Demux; Gaussian mixture with multi-sample-multiplet accounting.
- Klein 2023, Bioinformatics 39(8):btad481 - demuxmix; negative-binomial regression mixture robust to staining differences.
