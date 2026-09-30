# Enhancer-gene linking usage guide

## Prerequisites

```bash
# ABC (scripts/run_abc.sh calls its workflow/scripts/*.py); tested on tag v1.1.2 and main
git clone https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
# ABC's own workflow/envs/abcenv.yml did not solve with micromamba 2 (c-ares/libdeflate conflict, 2026-09-30).
# Tested package set instead:
micromamba create -n abc -c conda-forge -c bioconda python=3.10 "numpy<2" pandas pyranges hic-straw \
    macs2 samtools bedtools "snakemake=7" "pulp=2.7" pysam pybigwig scipy click
# macs2 2.2.9.1 fails to import on glibc >= 2.31 (undefined symbol __log_finite); either preload a
# shim for __log_finite or call peaks with macs3 (same options) and pass them as PEAKS.

# ENCODE-rE2G (runs ABC as a submodule)
git clone --recurse-submodules https://github.com/EngreitzLab/ENCODE_rE2G
# dependencies as in workflow/envs/encode_re2g.yml; unpinned `pip install snakemake` installs
# snakemake 9, which the pipeline does not support: keep snakemake>=7,<8 and pulp<2.8.
# Its pins that matter: bedtools=2.29.2, scikit-learn=1.2.1.

# HiC-Pro and FitHiChIP are not on bioconda: install HiC-Pro via its environment.yml
# (github.com/nservant/HiC-Pro) and FitHiChIP from github.com/ay-lab/FitHiChIP (Docker/Singularity)
```

Cicero and its R dependencies are installed by the atac-seq/co-accessibility Skill; `monocle3` is not a Bioconductor package, so `BiocManager::install('monocle3')` cannot resolve.

Inputs:
- Accessibility: DNase or ATAC BAM (ABC's docs also list ATAC tagAlign, not run here); replicates comma-separated
- H3K27ac ChIP-seq BAM matched to the accessibility sample (optional; selects the H3K27ac threshold row)
- Hi-C/Micro-C: `.hic` file or URL at 5 kb (`HIC_TYPE=hic`), or a directory for `avg`/`juicebox`/`bedpe`; none uses the powerlaw
- Gene annotation BED6 and blocklist (the ABC repository ships hg38 references and the K562 quantile-normalisation reference)

Example, K562 chr22 from ABC's `example_chr`, real ENCODE Hi-C read by range:

```bash
export ABC_REPO=$PWD/ABC-Enhancer-Gene-Prediction
ACCESS_TYPE=DHS OUTDIR=abc_out CELL_TYPE=K562 \
  ACCESS_BAM=$ABC_REPO/example_chr/chr22/ENCFF860XAE.chr22.sorted.se.bam \
  H3K27AC_BAM=$ABC_REPO/example_chr/chr22/ENCFF790GFL.chr22.sorted.se.bam \
  GENES_BED=$ABC_REPO/example_chr/chr22/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.hg38.bed \
  TSS_BED=$ABC_REPO/example_chr/chr22/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.hg38.TSS500bp.bed \
  HIC_TYPE=hic HIC_FILE=https://www.encodeproject.org/files/ENCFF621AIY/@@download/ENCFF621AIY.hic \
  bash scripts/run_abc.sh
```

This reproduces ABC's own chr22 result: 17,732 candidate regions, 1,674,535 expressed-gene pairs, 1,353 links at threshold 0.027.

## Example requests

- "Run ABC on K562 ATAC + H3K27ac + 5 kb Hi-C; use the calibrated threshold and drop promoter links."
- "Use ENCODE-rE2G on my B-cell DNase + H3K27ac with cell-type Hi-C; report links above the model threshold and cross-check against ABC."
- "I only have ATAC + H3K27ac. Run ABC with the powerlaw contact and note the degraded performance."
- "Run ABC and ENCODE-rE2G, intersect them, and flag HiChIP H3K27ac loop support; report the triple-overlap as high-confidence and the union as exploratory."
- "Compare my K562 ABC predictions against the ENCODE-rE2G CRISPR benchmark; report sensitivity and specificity at the calibrated threshold."
- "For a GWAS lead SNP in a non-coding region, find the overlapping ATAC + H3K27ac enhancer and its likely ABC target gene."

## Agent procedure

1. Verify inputs and choose the method by data availability (SKILL.md table).
2. ABC: `scripts/run_abc.sh` (peaks, candidate regions, Activity, Contact, ABC score, calibrated threshold). ENCODE-rE2G: write the biosamples TSV and let the pipeline pick the model. HiChIP: FDR < 0.05 and count >= 5, then intersect with ABC.
3. Cross-validate against CRISPR benchmark data where available.
4. Report intersections as high-confidence, unions as exploratory, and document Hi-C proxies.
