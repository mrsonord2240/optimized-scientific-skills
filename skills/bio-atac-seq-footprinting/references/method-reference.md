# TF footprinting method reference

Conditional detail for [`../SKILL.md`](../SKILL.md): tool choice, bias models, per-TF failure modes, interpretation, databases, errors, and citations.

## Tested-with versions

Executed on ENCODE GM12878 and K562 ATAC chr1 slices (hg38) and 10x PBMC 5k scATAC fragments (chr1:1-30 Mb) on 2026-09-30: TOBIAS 0.17.5, RGT HINT-ATAC 1.0.2, pyDNase (Wellington) 0.3.0, scPrinter 1.2.0 with tangermeme 0.4.4 and snapatac2 2.8.0, samtools 1.19.2, bedtools 2.31.1, pyBigWig 0.3.26, deepTools 3.5.5. The chr1 slices hold about 0.5-5M reads, so results verify that the tools run and separate bound from unbound sites, not statistical power at full depth. Check installed versions with `<tool> --version` and `--help` before relying on flags.

## Tool taxonomy

| Tool | Bias model | Scoring | Min depth | Strength | Fails when |
|------|-----------|---------|-----------|----------|------------|
| TOBIAS (BINDetect) | +/-12 bp k-mer window (`--k_flank` 12), dinucleotide weight matrix (DWM) | Two-step: continuous footprint score then motif-anchored bound/unbound classification | >= 50M nuclear reads | Mature, peer-reviewed (Bentsen 2020), differential support, modular pipeline | Below 50M reads; sequencing errors near motif inflate background |
| HINT-ATAC | Hidden-Markov + dinucleotide bias correction | HMM emits open/footprint/closed states; calls ranked footprints | >= 50M | Single-step; integrates motif matching; handles DNase too | Less control over individual stages; HMM occasionally over-segments |
| Wellington (pyDNase) | DNase-developed; `-A` ATAC mode | Cleavage-rate Poisson Z-score | >= 50M (DNase >= 80M) | Original footprinting framework; well-validated for DNase | Designed for DNase; ATAC-specific bias not corrected as carefully |
| PIQ | Bayesian latent variable on cut sites | Genome-wide PWM scan + cleavage profile | Not established here | Per-TF posterior probabilities | Outdated; not installed or run here |
| scPrinter | Pretrained Tn5 bias model plus multi-scale footprint scores | Footprint score at each scale (mode) around region centres | 50M (bulk); per-cluster: see the usage guide | Multi-scale; can resolve TF families with different footprint sizes | Newer tool; fragile install pins; GPU recommended; first import downloads models |

Command surfaces: `TOBIAS ATACorrect` -> `TOBIAS ScoreBigwig` (formerly `FootprintScores`) -> `TOBIAS BINDetect`; `rgt-hint footprinting --atac-seq` (HINT-ATAC, single step); `wellington_footprints.py -A` (DNase-derived, with an ATAC mode); Python `scprinter` (multi-scale; Hu 2025 Nature; bulk and per-cluster route in `scripts/scprinter_footprint.py`).

Methodology evolves; verify against the current Bentsen 2020, Karabacak Calviello 2019, and scPrinter (Hu 2025) benchmarks.

## Decision tree by goal

| Goal | Recommended pipeline |
|------|---------------------|
| Identify all TFs differentially bound between two conditions | TOBIAS ATACorrect (per condition) -> ScoreBigwig -> BINDetect with `--cond-names` |
| Find the strongest single-TF binding (e.g., CTCF) | TOBIAS PlotAggregate over JASPAR CTCF motif sites; verify the central dip |
| Per-cluster footprinting (scATAC) | `scripts/scprinter_footprint.py --groups`; tested depth and limits in the usage guide |
| Multi-scale TF activity (short and long simultaneously) | scPrinter (`scripts/scprinter_footprint.py`) |
| Differential nuclear-receptor binding | TOBIAS pooled-replicate footprints + ChIP cross-validation; ATAC alone often misses transient binding |
| Plant / non-model organism | TOBIAS or HINT-ATAC with custom motifs (for example CIS-BP); bias model retrained from genomic background |
| Single condition, find bound TFs | TOBIAS or HINT-ATAC |
| Lower depth (< 50M) | Pool replicates first; treat calls as exploratory |

## Tn5 bias and why correction matters

Tn5 inserts preferentially at certain k-mers (Karabacak Calviello 2019 measured the protocol-specific 6-mer insertion-bias model used for bias correction). The preference is reproducible and biologically uninteresting. Without correction, every "TF footprint" near a high-bias k-mer reads as occupancy, and regions with low-bias flanks but real binding may show no dip relative to corrected expectation. Symptom: aggregates at motifs of biased composition show a dip or peak that follows the bias track rather than a TF, so inspect the `_bias.bw` and `_expected.bw` profiles at the same sites. Fix: ATACorrect (TOBIAS), seqOutBias (Martins 2018), or HINT's dinucleotide model, which subtract the Tn5-expected per-base profile from observed cleavage; after correction the V-shape should persist only at TF-bound sites.

### Bias correction alternatives

| Method | Approach | When to use |
|--------|---------|-------------|
| TOBIAS ATACorrect | +/-12 bp k-mer window, dinucleotide weight matrix (DWM) | Default for most ATAC; fast |
| chromBPNet bias model (Pampari 2024) | CNN trained on naked-DNA control or k-mer baseline | Best when sequence context is complex; handles low-complexity flanks |
| seqOutBias (Martins 2018) | Genome-wide k-mer frequency scaling (observed vs expected cut counts) | Independent of footprinting tool; works upstream |
| HINT-ATAC dinucleotide | HMM-integrated dinucleotide bias | Built into HINT pipeline; less control |
| Naked-DNA empirical | Sequence Tn5 on protein-free DNA | Gold standard for non-model organisms; expensive wet-lab |

For non-model organisms (no published Tn5 bias model), a naked-DNA control is preferred. chromBPNet's bias model (Pampari 2024) is reported to handle low-complexity sequence contexts better than TOBIAS; it was not run here, see the deep-learning ATAC Skill.

## In silico variant effect at footprinted motifs

When a GWAS-fine-mapped or rare variant lies inside a TOBIAS-bound motif site, sequence-based deep-learning models (chromBPNet, Enformer) predict per-base accessibility for ref versus alt; combined with footprint evidence this gives a mechanistic hypothesis ("variant disrupts binding of TF X at enhancer Y"). Run BINDetect to identify bound sites, score variants in bound sites with chromBPNet for ref/alt log2FC (|log2FC| > 1 supports functional disruption), and cross-reference allele-specific accessibility for observed evidence.

## Per-TF footprinting failure modes

The same tool can give a clean V-shape for one TF family and noise for another.

- **CTCF, the usual positive control.** Stable binder with a long, deep footprint, strong sequence specificity, and a clean dip with flanking cleavage shoulders. In the chr1 GM12878 test the corrected CTCF aggregate showed a central dip with flanking peaks at 513 bound sites and a flat profile at 690 unbound sites (1203 sites in total); the uncorrected signal at the same bound sites also dipped, so the dip alone does not validate the correction. Over all 1203 sites, the aggregate of ATACorrect's bias-only expectation (`_expected.bw`) correlated with the uncorrected aggregate (Pearson r 0.65 in GM12878, 0.40 in K562) but not with the corrected one (-0.14, -0.22). With the corrected track replaced by the uncorrected one, the bound-site dip stayed as deep while r rose to 0.65 and 0.40, so the `run_tobias.sh` check (fail above 0.2) exited 4. Always check CTCF first; if it is shallow, bias correction or depth is the likely problem, not the biology.
- **Nuclear receptors (ER, AR, GR), transient binding.** Residence time is short compared with CTCF, so the average ATAC sample captures a low binding probability per allele. The aggregate footprint is shallow or absent despite ChIP-seq peaks. Use scprinter's multi-scale model, limit to ChIP-validated sites, or pool replicates. Do not interpret absence of footprint as absence of binding.
- **Pioneer TFs (FOXA1, GATA, OCT4), half-site footprint.** They bind one DNA face while the other is on the histone octamer, giving an asymmetric V with one taller shoulder. Inspect aggregates per motif strand; no strand-specific calling option was found in rgt-hint 1.0.2. Treat asymmetry as possibly biological, not automatically artefactual.
- **AP-1 (FOS, JUN), composite footprint.** JASPAR entries are degenerate across heterodimer combinations (FOS+JUN, FOS+JUNB, JUNB+JUNB), so scoring averages over them. Use specific HOCOMOCO motifs per heterodimer when distinguishing matters; otherwise accept the composite call.
- **ZBTB / BTB-zinc finger, dynamic or unfootprintable.** Dynamic kinetics and cofactor-mediated stabilization make steady-state occupancy variable; some (ZBTB16, BCL6) do not yield reliable footprints despite genuine ChIP-seq binding. Document the failure and use ChIP-seq.
- **Forkhead / homeobox (FOX, HOX), short footprints under 8 bp.** Footprint extent matches motif length, at the resolution limit of Tn5 (about 4 bp positional uncertainty). Use multi-scale scoring (scprinter), aggregate over thousands of sites, and do not rely on per-site calls.

ATAC footprints do not replicate ChIP-seq for all TFs: concordance is better for stable binders such as CTCF and poorer for transient ones such as nuclear receptors (no numeric expectation is established here).

## Reading BINDetect differentials

| BINDetect output | Interpretation |
|------------------|----------------|
| `cond1_cond2_change` > 0, low pvalue | TF more bound in cond1 |
| `cond1_cond2_change` < 0, low pvalue | TF more bound in cond2 |
| Both `cond1_bound` and `cond2_bound` near 0 | Motif present but no footprint either condition; TF likely not active |
| `cond1_bound` >> `cond2_bound` but change small | High dynamic range; differential per-site rather than aggregate |

Output columns: `output_prefix`, motif info, condition counts (`cond1_bound`, `cond2_bound`), `cond1_mean_score`, `cond2_mean_score`, `cond1_cond2_change`, `cond1_cond2_pvalue` (one-sided per direction). Motif variants (for example the CTCF IDs MA0139.2, MA1929.2, MA1930.2) are separate rows and directories. Observed in the chr1 GM12878 (cond1) versus K562 (cond2) test: IRF4 +0.59, GATA1 -0.39, CEBPA +0.09, EBF1 +0.06; values above 1.0 would warrant investigation.

## Reconciling TOBIAS and HINT-ATAC

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| TOBIAS calls binding, HINT does not | TOBIAS more sensitive; HINT's HMM filters edge calls | Trust if motif is canonical; suspect for novel/weak motifs |
| HINT calls binding, TOBIAS does not | HINT's HMM occasionally over-segments and reports spurious calls | Verify by aggregate footprint at the called sites |
| Both call same TF as differential but opposite directions | Different bias correction model; different bound/unbound thresholds | Re-check ATACorrect output; one bias model may be miscalibrated |
| Both flat | Library too shallow; chromatin too closed at motif sites | Pool replicates; consider scprinter multi-scale |

Concordance measure (`scripts/site_concordance.sh`): take TOBIAS bound and unbound motif sites inside the regions the second method covered, and report the fraction of each set that overlaps the second call set (HINT-ATAC or Wellington footprints, or ChIP-seq peaks), plus the reverse fraction (second-set calls overlapping a bound site). Higher confidence means the bound-site fraction is clearly above the unbound-site fraction; no universal percentage cutoff is established, and a fixed 50% would have rejected the real results below. Report single-tool calls as exploratory.

Measured on GM12878 chr1:10-13 Mb (60 peaks, about 0.5M fragments): HINT-ATAC footprints overlapped 58/143 bound sites (41%) and 0/31 unbound sites, with 15% of its 358 footprints on a bound site; Wellington `-A` gave 12% and 0/31 (29% of its footprints), Wellington without `-A` 6% and 0/31. Low depth and site-level overlap make these fractions a check of agreement, not a benchmark.

## Motif database choice

| Database | Coverage | Format | Notes |
|----------|---------|--------|-------|
| JASPAR 2024 CORE vertebrates | ~880 vertebrate non-redundant motifs (curated, experimentally derived) | JASPAR PFM, MEME, etc. | Default for vertebrate ATAC |
| HOCOMOCO v12 | ~1443 curated motifs (v12 CORE) | JASPAR PFM | Best for resolving paralogues; secondary motif subtypes per TF |
| CIS-BP 2.0 | ~80,000 motifs across 1000+ species | PWM, .meme | Broadest coverage including non-model species |
| MEME-CHIP / homer | Custom from peaks | .meme, .motif | When de novo motif needed |

JASPAR is conservatively curated; HOCOMOCO is comprehensive for human/mouse with per-motif quality scores (A/B/C/D); CIS-BP excels for non-model organisms.

## Common errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Aggregate footprint inverted (peak instead of dip) | No bias correction; or wrong genome FASTA | Run ATACorrect; verify FASTA matches BAM build |
| BINDetect reports zero bound sites | Default cutoff too stringent; or peakset too narrow | Raise `--bound-pvalue` (default 0.001) and inspect; verify peaks include where binding is expected |
| TOBIAS ATACorrect out of memory | Genome FASTA huge or many cores | Reduce `--cores`; use `samtools faidx` to confirm the FASTA index exists |
| Differential score noisy / random | Per-condition bias correction inconsistent | Re-run ATACorrect with identical peakset and blacklist for each |
| Empty motif file warning | JASPAR PFM format mismatch | Use MEME suite to convert; TOBIAS expects JASPAR format |
| HINT-ATAC reports many tiny footprints | Default HMM over-segments | Use `--organism` explicitly; make the regions BED the consensus peakset, not raw peaks |
| `rgt-hint` fails with `FileNotFoundError ... data.config` | The bioconda `rgt` package ships no data | Set up RGT data (usage guide) and export `RGTDATA` |
| Wellington run without ATAC mode | `-A` omitted | Use `wellington_footprints.py -A`; paired-end ATAC BAMs run to completion in pyDNase 0.3.0 |
| Per-site footprint is V-shape but aggregate is flat | Mixing strands; some motifs on - strand | Aggregate function should handle strand; verify input motif strand column |

## References

- Buenrostro JD et al 2013 Nat Methods 10:1213 (ATAC-seq protocol)
- Karabacak Calviello A et al 2019 Genome Biol 20:42 (protocol-specific Tn5/DNase bias modeling for footprinting)
- Bentsen M et al 2020 Nat Commun 11:4267 (TOBIAS framework, benchmark)
- Li Z et al 2019 Genome Biol 20:45 (HINT-ATAC)
- Piper J et al 2013 NAR 41:e201 (Wellington / pyDNase)
- Sherwood RI et al 2014 Nat Biotechnol 32:171 (PIQ)
- Hu Y et al 2025 Nature 638:779 (scPrinter/PRINT; multi-scale single-cell footprinting)
- Martins AL et al 2018 NAR 46:e9 (seqOutBias)
- Castro-Mondragon JA et al 2022 NAR 50:D165 (JASPAR 2022, recently 2024)
- Vorontsov IE et al 2024 NAR 52:D154 (HOCOMOCO v12)

## Related bioSkills topics (not bundled)

atac-seq peak calling (input peakset), atac-seq QC (confirm depth >= 50M), single-cell ATAC (scprinter), motif deviation (complementary chromVAR), deep-learning ATAC (chromBPNet bias model and variant effects), allele-specific accessibility, chip-seq peak annotation (cross-validation), sequence-manipulation motif search, gene-regulatory-networks SCENIC regulons.
