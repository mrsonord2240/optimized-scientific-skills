# Nucleosome positioning method reference

## Nucleosome physics for ATAC

A nucleosome wraps about 147 bp DNA in 1.65 turns. Adjacent nucleosomes are separated by 20-50 bp linker; the mean **nucleosome repeat length (NRL)** is species-dependent:

| Cell type / organism | NRL | Notes |
|----------------------|-----|-------|
| Yeast S. cerevisiae | 165 bp | Tightly packed; less linker |
| Drosophila S2 | 175-185 bp | |
| Mouse ES cells | 188-196 bp | |
| Human HEK293 / K562 | 196-200 bp | Standard somatic |
| Human cortical neurons | 211 bp | Longer linker |
| Sperm chromatin | 240-250 bp | Tight packaging via protamines |
| Active gene bodies | -10 bp shorter than genome avg | Active transcription disrupts |

NRL determines fragment-size peak positions. The ATAC mono-nucleosome peak is at NRL, not 147 bp (147 bp is the protected length; ATAC fragments span nucleosome plus linker). The di-nucleosome peak is at 2x NRL minus a small overlap.

## Fragment-size classes (Buenrostro 2013, refined)

| Class | Fragment range | Origin | Use |
|-------|---------------|--------|-----|
| Sub-nucleosomal / NFR | < 100 bp | Two Tn5 cuts in naked accessible DNA | TF binding, footprinting |
| Mono-nucleosomal | 180-247 bp | Tn5 cuts on each side of one nucleosome | Nucleosome positioning |
| Di-nucleosomal | 315-473 bp | Tn5 cuts span two nucleosomes | Phasing, NRL estimation |
| Tri-nucleosomal | 558-615 bp | Three nucleosomes | Heterochromatin / phasing |
| > 700 bp | Rare | Often artefact (chimeric); discard | -- |

The 180-247 bp mono window is the Buenrostro 2013 convention; ATACseqQC uses 180-250. Adjust for the target organism's NRL.

## V-plot interpretation

V-plots (fragment size vs position) are diagnostic. The x-axis is position relative to a feature (TSS, motif center); the y-axis is fragment size.

| Pattern | Visual | Meaning |
|---------|--------|---------|
| V (apex at center, low size at center, increasing flanks) | Classic V | TF or NFR at center, flanking nucleosomes |
| W (two V's flanking center) | W-shape | NFR at center plus +1 / -1 nucleosomes |
| Inverted V (peak at center) | Mountain | Fragment fully enclosed at feature; e.g. nucleosome-bound TF |
| Flat band at 200 bp | Horizontal line | Constitutive nucleosome (no positioning relative to feature) |
| 10.4 bp helical phasing on V | Sub-peaks at 50, 60, 70, 80 bp size | Tn5 helical preference visible; high-quality library |

V-plots are the primary diagnostic for whether positioning analysis will succeed. Flat bands mean no positioning information; classic V/W patterns mean positioning is recoverable. Use the script `scripts/vplot.py`.

## Algorithmic taxonomy

| Tool | Method | Resolution | Strength | Fails when |
|------|--------|-----------|----------|------------|
| NucleoATAC | Cross-correlation with idealized V-plot template; per-base occupancy + nucleosome calls | Single-bp | ATAC-specific; provides occupancy + fuzziness | Unmaintained (last release 0.3.4, last commit 2019-03); Python 2.7 only; struggles on chromatin without clear NRL |
| ATACseqQC | Fragment-size split + Tn5-shifted GAlignments + V-plot from BAM | Region-level | R/Bioconductor; integrates with TxDb / motif analysis | No per-base nucleosome calls; visualization-focused |
| DANPOS3 | Smoothing + peak call on cleavage signal; tested on MNase, ATAC, DNase | ~50 bp | Differential mode (`dpos`); MNase legacy | Designed for MNase-Seq; ATAC adaptation needs careful parameter tuning |
| scprinter | CNN multi-scale; resolves co-occurring TF + nucleosome footprints | Single-bp | Modern; single-cell aware; multi-scale | Newer; benchmarks evolving; GPU recommended |
| custom (pysam V-plot) | Fragment counting + 2D density | Region-level | Maximally flexible; reproducible | Requires manual calling logic; slow |

scprinter is named as an alternative but this Skill ships no tested scprinter workflow; follow its own documentation. Verify against Schep 2015 (NucleoATAC), Chen 2013 (DANPOS), Hu 2025 (scPrinter) before locking pipelines.

## +1 nucleosome calling

The +1 nucleosome (first nucleosome downstream of the TSS, bordering the NFR) is the most-studied positioning feature; its position relative to TSS determines transcription initiation kinetics.

Canonical +1 position: +50 to +60 bp from TSS in metazoa (yeast and other organisms differ; no single figure is asserted here, see Mavrich 2008); varies by gene type (Pol II vs Pol III, housekeeping vs developmental).

```bash
# 1. Define gene-body intervals
bedtools slop -i genes.bed -g chrom.sizes -l 200 -r 1000 > gene_bodies.bed

# 2. Run NucleoATAC
nucleoatac run --bed gene_bodies.bed --bam sample.dedup.bam --fasta genome.fa \
    --out tss_nuc/ --cores 8

# 3. The first nucleosome downstream of each TSS in nucpos.bed is +1
#    (only after the nucpos check in the NucleoATAC section below)
```

Failure to detect a clear +1 peak in the aggregate V-plot suggests the TSS annotation is wrong or the library is over-transposed.

## Per-tool failure modes

### NucleoATAC: region size and depth dependence

- Trigger: short region BED (< 1 kb per region); shallow library (< 25M nuclear reads).
- Mechanism: NucleoATAC fits an idealized V-plot template per region; short regions give too few fragments for stable correlation and shallow data gives noisy templates.
- Symptom: no nucleosome calls in shallow regions; occupancy track flat at zero.
- Fix: use regions >= 500 bp; merge adjacent peaks via bedtools; require >= 30M nuclear reads.

### NucleoATAC: maintenance status and empty nucpos

- Trigger: installing NucleoATAC now, or `nucleoatac run` exits 0 but `out.nucpos.bed.gz` is empty.
- Mechanism: 0.3.4 is the only release (last commit 2019-03) and installs only on Python 2.7; `pip install nucleoatac` on Python 3 fails with "Python version must be 2.7!". Its compiled `multinomial_cov.calculateCov` declares `cdef DTYPE_t value` without initialising it, so z-scores are NaN or garbage and every candidate is rejected (`--min_z 0` cannot help); occupancy, NFR and occupancy-peak outputs are unaffected.
- Fix: install per the usage guide (bioconda Python 2.7 build, initialise `value = 0` in the installed `multinomial_cov.pyx`, rebuild with `cythonize -i`). This edits third-party code inside your NucleoATAC environment; record it in the methods. Always check after a run that `out.nucpos.bed.gz` is non-empty; an empty file on regions of 500 bp or more means the defect, not absent nucleosomes. If you cannot patch, report nucleosome calls as unavailable and use occupancy/NFR outputs or DANPOS3. Consider DANPOS3 for new projects.
- Tested outcome: GM12878 rep1+rep2 (0.9M pairs, chr1:1-30 Mb), 209 TSS-flank regions: 178 nucpos calls, 6 redundant; median first call downstream of a TSS was +121 bp (IQR 89-225) at this shallow depth.

### ATACseqQC factorFootprints: asymmetric nucleosome flanks

- Trigger: pioneer-factor binding sites where one face is on a nucleosome.
- Mechanism: factorFootprints assumes symmetric flanking nucleosomes; pioneer TFs (FOXA1, GATA) have a nucleosome on one side only, giving asymmetric output.
- Symptom: single shoulder in flanking signal; unbalanced V-plot.
- Fix: treat asymmetry as biological signal, not artefact. For pioneers, use stranded analysis.

### DANPOS dpos with default parameters

- Trigger: running `danpos dpos` with MNase defaults on ATAC.
- Mechanism: DANPOS3 smoothing window and peak-calling defaults are tuned for MNase signal (smoother coverage); ATAC's sharper signal needs `--smooth_width 80 --width 145` or similar, otherwise calls are over-smoothed.
- Fix: use ATAC-tuned parameters (see below and DANPOS docs), or use NucleoATAC instead.

### Mono-nucleosome filter window mis-set

- Trigger: strict 147 bp filter for the mono-nuc fraction; using 100-180 bp instead of 180-247.
- Mechanism: mono-nuc fragments are 180-247 bp because they span the nucleosome and a linker; tighter filtering excludes real signal.
- Symptom: mono-nuc count much lower than expected (< 30% of NFR count).
- Fix: use Buenrostro 2013 windows: NFR < 100, mono 180-247, di 315-473.

## Decision tree by goal

| Goal | Recommended workflow |
|------|---------------------|
| Per-base nucleosome occupancy track | NucleoATAC (with caveat about maintenance); or scprinter |
| V-plot at TSS or motif center | ATACseqQC vPlot (motif-centered; needs a pfm and binding sites) or `scripts/vplot.py` |
| Differential nucleosome positioning between conditions | DANPOS3 dpos |
| +1 nucleosome calling at all genes | NucleoATAC + post-process to first nuc downstream of TSS |
| Single-cell nucleosome positioning | scprinter |
| Quick fragment-size QC plot | ATACseqQC fragSizeDist |
| NRL estimation | `scripts/estimate_nrl.py` (approximate); autocorrelation or Fourier on cumulative cleavage coverage for precision |
| Nucleosome-aware peak calling | MACS3 hmmratac (peak-calling skill) |

## Differential positioning (DANPOS3 dpos)

```bash
# Compare control vs treatment nucleosome positions.
# The sample pair is the POSITIONAL argument (a:b means a minus b); -b is for background/input to
# subtract, and -c specifies a read-count to normalize to (an integer, NOT a control BAM path).
danpos dpos condition2.bam:condition1.bam \
    -o danpos_diff/ \
    --paired 1 \
    --smooth_width 80
```

The differential table is `condition2-condition1.positions.integrative.xls`. Its shift column is `treat2control_dis` (signed bp; `diff_smt_loca` is a coordinate) and the FDR that tracks position shifts is `point_diff_FDR`. `smt_diff_FDR` tests summit height and removed every planted shift in testing. The table has no event-type label column; classify gained, lost, shifted, and fuzziness-changed nucleosomes from the value and FDR columns yourself. ENCODE has no official threshold; require |shift| >= 30 bp and `point_diff_FDR` < 0.05:

```bash
awk -F'	' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next}
     $c["treat2control_dis"]!="NA" && $c["point_diff_FDR"]!="NA" &&
     ($c["treat2control_dis"]>=30 || $c["treat2control_dis"]<=-30) && $c["point_diff_FDR"]<0.05'     danpos_diff/condition2-condition1.positions.integrative.xls > shifted.tsv
```

Full ATAC-tuned DANPOS3 recipe:

```bash
# --width 145: summit-scan window (DANPOS -jw/--width; default 40)
# --smooth_width 80: smoothing kernel width (DANPOS -z; default 20, widened for ATAC)
# -jd 145: min distance between adjacent nuc calls (single-dash short flag)
# --pheight 1e-5: occupancy P-value cutoff (DANPOS -p; dpos default 0). -q/--height is the separate density cutoff
# --frsz 200: fragment size used (mono-nuc)
danpos dpos sample.bam \
    --paired 1 \
    --width 145 \
    --smooth_width 80 \
    -jd 145 \
    --pheight 1e-5 \
    --frsz 200 \
    --out danpos_out/
```

Verify exact flags with `danpos dpos --help`. The tested route is the conda `danpos3` package (3.2.4), which installs an executable named `danpos` and prints a 3.1.1 banner; a GitHub checkout (github.com/sklasfeld/DANPOS3) runs as `python danpos.py` with the same flags. The bioconda `danpos` package is DANPOS2.

Adapted from DANPOS3 docs for ATAC; `--smooth_width 80` widens the smoothing kernel to match ATAC's sharper signal versus MNase's broader cleavage. `-jd 145` (alternative `--distance 145`) enforces nucleosome spacing >= 145 bp (one nucleosome footprint).

## Histone variant detection from fragment size

- Trigger: suspected H2A.Z- or H3.3-containing nucleosomes; differential nucleosome composition between conditions.
- Mechanism: H3.3-H2A.Z double-variant nucleosomes are destabilized at active promoters (Jin 2009 Nat Genet 41:941-945). That paper does not report an ATAC fragment-size shift, so any size difference between H2A.Z-positive and -negative regions is an untested hypothesis to be measured, not assumed.
- Detection: compare the aggregate fragment-size distribution at H2A.Z ChIP-seq peaks versus H3K4me3-only peaks. ATAC alone cannot call H2A.Z; H2A.Z ChIP-seq is needed for ground truth.

```python
import numpy as np

# Per-region mean fragment size (bam is an open pysam.AlignmentFile; region = (chrom, start, end))
def region_frag_size(bam, region):
    sizes = [abs(r.template_length) for r in bam.fetch(*region)
             if r.is_proper_pair and r.is_read1 and 100 < abs(r.template_length) < 300]
    return np.mean(sizes) if sizes else np.nan

# Compare H2A.Z-positive vs H2A.Z-negative TSSs
```

## Long-read single-molecule chromatin

| Method | Tech | Resolution | Strength |
|--------|------|-----------|----------|
| Fiber-seq (Stergachis 2020) | PacBio HiFi + DNA methylation footprinting | Per-molecule single-bp | Reads continuous chromatin fiber up to 20 kb; resolves haplotype-specific positioning |
| NanoNOMe (Lee 2020 Nat Methods 17:1191-1199) | Nanopore + GpC methyltransferase | Per-molecule single-bp | Same single-molecule resolution, cheaper than PacBio |

Fiber-seq detects nucleosome occupancy directly per chromatin molecule (no aggregation) and resolves cell-cycle-dependent and stochastic positioning that bulk ATAC averages out. Prefer it for fine-structure analysis of regulatory elements. For most labs short-read ATAC plus NucleoATAC remains primary.

## Nucleosome fuzziness

Fuzziness measures how sharply positioned a nucleosome is across cells, defined as the standard deviation of per-cell nucleosome center positions.

| Fuzziness range | Interpretation |
|-----------------|----------------|
| < 20 bp | Sharply positioned (rare in metazoa; common at +1 in yeast) |
| 20-50 bp | Standard well-positioned |
| 50-100 bp | Fuzzy; constitutive but non-stable |
| > 100 bp | Effectively unpositioned |

These are field-convention bands (NucleoATAC / DANPOS practice); no single primary paper prescribes them, so verify against tool documentation when reporting. NucleoATAC reports per-nucleosome fuzziness in column 13 of `.nucpos.bed`.

## Common errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `pip install nucleoatac` says Python version must be 2.7 | Python 3 environment | Python 2.7 conda env (usage guide) |
| Empty .nucpos.bed output with regions >= 1 kb | Uninitialised `calculateCov` accumulator (NaN z-scores) | Apply the `value = 0` rebuild (section above) |
| Empty .nucpos.bed output on tiny regions or shallow data | Region BED too short or library too shallow | Verify region size >= 500 bp; depth >= 30M |
| V-plot shows horizontal band, no V | No positioning info; library over-transposed or wrong feature center | Check feature BED; verify TSS positions |
| Mono-nuc count very low | Wrong fragment-size window (100-180 instead of 180-247) | Use Buenrostro windows |
| factorFootprints asymmetric | Pioneer TF; biological | Treat as signal, not artefact |
| DANPOS calls many shifts | MNase parameters used on ATAC | Tune `--smooth_width 80 --width 145` for ATAC |
| splitGAlignmentsByCut error in ATACseqQC | BAM is single-end | Mono-nuc analysis requires paired-end |
| +1 nucleosome not visible at TSS aggregate | TSS list mixes coding and non-coding strands; or wrong genome build | Restrict to protein-coding TSSs in matched build |

## References

- Schep AN et al 2015 Genome Res 25:1757 (NucleoATAC)
- Chen K et al 2013 Genome Res 23:341 (DANPOS)
- Buenrostro JD et al 2013 Nat Methods 10:1213 (ATAC fragment-size classes)
- Ou J et al 2018 BMC Genomics 19:169 (ATACseqQC)
- Hu Y et al 2025 Nature 638:779 (scPrinter/PRINT; multiscale footprints)
- Mavrich TN et al 2008 Nature 453:358 (+1 nucleosome positioning)
- Voong LN et al 2016 Cell 167:1555-1570 (high-resolution chemical nucleosome mapping)
- Jin C et al 2009 Nat Genet 41:941 (H3.3/H2A.Z double-variant nucleosome instability)
- Teif VB et al 2012 Nat Struct Mol Biol 19:1185 (NRL variation across cell types)

## Related Skills (bioSkills paths)

- atac-seq/atac-qc, atac-seq/atac-peak-calling, atac-seq/footprinting, atac-seq/single-cell-atac
- chip-seq/peak-annotation, alignment-files/bam-statistics
