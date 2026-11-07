# Interpretation



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
