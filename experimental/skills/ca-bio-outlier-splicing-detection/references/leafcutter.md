# LeafcutterMD details


**Goal:** Detect outlier intron usage relative to a control panel without annotation dependence.

**Approach:** Run LeafcutterMD (LeafCutter's Dirichlet-multinomial outlier mode for Mendelian disease) against the control panel. `leafcutter_cluster_regtools.py` is in the leafcutter repository (`clustering/`), `leafcutterMD.R` in `scripts/`; neither is installed by `BiocManager`. BAMs must be indexed and carry XS strand tags (regtools writes strand `?` otherwise and the clustering step drops those junctions).

```bash
for bam in *.bam; do
    regtools junctions extract -a 8 -m 50 -s XS "$bam" -o "${bam%.bam}.junc"
done

ls *.junc > juncfiles.txt
python /path/to/leafcutter/clustering/leafcutter_cluster_regtools.py -j juncfiles.txt -o leafcutter -m 50 -l 500000

/path/to/leafcutter/scripts/leafcutterMD.R \
    --num_threads 4 \
    --output_prefix patient_outlier \
    leafcutter_perind_numers.counts.gz
```

Outputs: `patient_outlier_clusterPvals.txt` (cluster x sample), `_pVals.txt` (intron x sample; row names `chr:start:end:clu_N_strand` give the coordinates of a cluster), `_effSize.txt` (intron x sample). The p-values are **raw**; control the FDR yourself:

```r
p <- read.table('patient_outlier_clusterPvals.txt', header = TRUE, sep = '\t', check.names = FALSE)
long <- data.frame(cluster = rep(rownames(p), ncol(p)), sampleID = rep(colnames(p), each = nrow(p)),
                   p = unlist(p, use.names = FALSE))
long <- long[!is.na(long$p), ]
long$q <- p.adjust(long$p, method = 'BH')
hits <- long[long$q < 0.05, ]
```

LeafcutterMD (Jenkinson 2020 *Bioinformatics*) is junction-only: intron retention and pseudoexon events produced no cluster p-value in testing. Useful when FRASER's model fits poorly or when novel-junction sensitivity matters.
