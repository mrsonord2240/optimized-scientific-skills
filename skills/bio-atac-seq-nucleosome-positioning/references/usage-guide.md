# Nucleosome positioning usage guide

## Prerequisites

```bash
# NucleoATAC 0.3.4 installs only on Python 2.7 (pip on Python 3 fails with "Python version must be 2.7!").
# Use an isolated conda env, then apply the one-line fix for its NaN-z-score defect (see the method reference).
micromamba create -n nucleoatac -c conda-forge -c bioconda --override-channels python=2.7 nucleoatac "cython<3"
P=$(micromamba run -n nucleoatac python -c "import nucleoatac,os;print(os.path.dirname(nucleoatac.__file__))")
sed -i 's/cdef DTYPE_t value$/cdef DTYPE_t value = 0/' $P/multinomial_cov.pyx
(cd $P && CFLAGS="-I$(micromamba run -n nucleoatac python -c 'import numpy;print(numpy.get_include())')"     micromamba run -n nucleoatac cythonize -i multinomial_cov.pyx)

conda install -c bioconda samtools bedtools
# DANPOS3: `conda install -c bioconda danpos3` gives the `danpos` executable (bioconda `danpos` is DANPOS2);
# alternatively clone github.com/sklasfeld/DANPOS3 and run `python danpos.py`.
# scPrinter installs from source: git clone https://github.com/buenrostrolab/scPrinter && cd scPrinter && pip install ./
pip install pysam pyBigWig matplotlib scipy
```

```r
BiocManager::install(c('ATACseqQC', 'ChIPpeakAnno', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
                       'BSgenome.Hsapiens.UCSC.hg38', 'rtracklayer', 'Rsamtools', 'MotifDb', 'motifmatchr'))
```

Tested with R 4.4.3 / Bioconductor 3.20, NucleoATAC 0.3.4 (Python 2.7.15, numpy 1.16.5, cython 0.29.15), danpos3 3.2.4, pysam 0.24.1, samtools 1.24, pyBigWig 0.3.26. `nucleosome_analysis.R` takes about 6 minutes on 0.5M read pairs with a protein-coding TSS BED (13 minutes with all knownGene transcripts); on a chromosome-slice BAM most TSS rows in the heatmap are empty.

Inputs: deduplicated, MAPQ-filtered, chrM-stripped paired-end BAM with >= 30M nuclear reads.

## Example requests

- Per-region calling: "Run NucleoATAC on consensus peaks merged with bedtools (regions >= 500 bp) in the Python 2.7 NucleoATAC env; report nucpos.bed and the occupancy bigWig."
- V-plot: "Generate a V-plot at TSSs (+/- 1 kb) to assess whether positioning is recoverable."
- +1 nucleosomes: "Define gene intervals around each TSS covering the +1 nucleosome (about +50 to +60 bp in metazoa), run NucleoATAC, and report the first nucleosome downstream of each TSS."
- Differential: "Compare two conditions with DANPOS3 dpos using `--width 145 --smooth_width 80`; keep |treat2control_dis| >= 30 bp at point_diff_FDR < 0.05."
- NRL: "Estimate the NRL from fragment sizes and compare with the expected value for the cell type."
- Single-cell: "Run scprinter on the scATAC fragment file; aggregate nucleosome tracks by cluster."

## Agent sequence

1. Verify paired-end BAM.
2. V-plot at TSSs (or another feature); abort with a diagnosis if it is a flat band.
3. Define regions (merged consensus peaks, TSS-flanking, or genome-wide).
4. Choose a tool: NucleoATAC for ATAC-specific calls, DANPOS3 for differential, scprinter for single-cell.
5. Report positions, occupancy, and fuzziness; annotate +1/-1/+2 nucleosomes relative to TSS; optionally estimate NRL.

## Tips

- Always plot the V-plot first; a flat band makes downstream calling invalid.
- The mono-nucleosome window is 180-247 bp, not 147 bp.
- NRL is species- and cell-type-specific; verify the mono-nuc window against the fragment-size distribution.
- The +1 nucleosome is about 50-60 bp downstream of the TSS in metazoa (other organisms differ); its absence in the aggregate suggests wrong TSS annotation.
- DANPOS3 MNase defaults are wrong for ATAC.
- scprinter is GPU-recommended; CPU runs are slow on large datasets.
- Shifts under 30 bp are within fragment-size noise.
- Fuzziness 20-50 bp is well-positioned in metazoa; over 100 bp is effectively unpositioned.
