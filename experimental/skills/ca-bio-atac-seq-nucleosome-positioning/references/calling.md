# Calling

The source slop example below does not specify strand-aware behavior; do not treat its left/right bounds as upstream/downstream on both strands without adjusting them. Verify the input BED strand and coordinate semantics before +1 assignment.



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
