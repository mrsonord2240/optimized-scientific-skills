## leafcutter Differential Intron Usage

**Goal:** Detect differential intron-cluster usage annotation-free, capturing novel junctions and complex multi-junction events.

**Approach:** Extract junctions with regtools, cluster introns by shared splice sites, run cluster-level Dirichlet-multinomial test.

```bash
# Strand comes from the XS tag: the BAMs need it (STAR --outSAMstrandField intronMotif); without XS regtools writes strand '?' and leafcutter drops every junction
for bam in *.bam; do
    regtools junctions extract -a 8 -m 50 -s XS "$bam" -o "${bam%.bam}.junc"
done
ls *.junc > juncfiles.txt

# LEAFCUTTER = a clone of github.com/davidaknowles/leafcutter (the R package does not install these scripts on PATH)
python $LEAFCUTTER/clustering/leafcutter_cluster_regtools.py \
    -j juncfiles.txt \
    -o leafcutter \
    -m 50 \
    -l 500000
```

- `-m 50` (min reads per cluster) yields 0 clusters on shallow data (a ~100k-read chrX subset needed `-m 10`).
- Contigs other than `chr1..22,X,Y` / `1..22,X,Y` (e.g. a custom or scaffold reference) are dropped silently and give 0 introns: add `-k True`.
- The counts-table column names are the `.junc` file names without `.junc`.

```bash
# groups.txt: two columns, no header; column 1 must equal the counts-table column names (the .junc basenames, e.g. G1_rep1), column 2 the group
# Extra columns = confounders (see Confounder Handling)
Rscript $LEAFCUTTER/scripts/leafcutter_ds.R --num_threads 4 -i 3 -g 3 -c 10 \
    --exon_file gencode_exons.txt.gz \
    leafcutter_perind_numers.counts.gz groups.txt -o ds_results
```

```r
cluster_sig <- read.table('ds_results_cluster_significance.txt', header = TRUE, sep = '\t')
intron_effects <- read.table('ds_results_effect_sizes.txt', header = TRUE, sep = '\t')

sig_clusters <- subset(cluster_sig, p.adjust < 0.05)
```

**Sample-size flags.** The defaults `-i 5 -g 3 -c 20` stop with `The number of samples in the smallest group is less than min_samples_per_intron` whenever a group has fewer than 5 samples (they stopped at 2v2, 3v3 and 4v4; 5v5 ran). Set `-i` and `-g` to the smaller group size and `-c` (reads per sample) to about 10:

| Design | Flags | Planted 24 strong / 6 weak events found | Null (A vs A') |
|--------|-------|------------------------------------------|----------------|
| 2v2 | `-i 2 -g 2 -c 10` | 21 / 1 | 0 calls |
| 3v3 | `-i 3 -g 3 -c 10` | 22 / 3 | 0 calls |
| 4v4 | `-i 4 -g 4 -c 10` | 23 / 2 | 0 calls |
| 5v5 | `-i 5 -g 5 -c 10` | 24 / 4 | 0 calls |

Simulated data, p.adjust < 0.05 and max |ΔPSI| > 0.10; 4 of 26 calls at 2v2 and 5 of 30 at 3v3 were false, 1 of 26 at 4v4 and 1 of 29 at 5v5, so treat n<4 calls as needing a second tool. The `--min_samples_per_intron 5` "pre-filter" is the default that stops small designs; do not add it.

`--exon_file` (columns chr, start, end, strand, gene_name) only labels clusters with genes; build it with `Rscript $LEAFCUTTER/scripts/gtf_to_exons.R annotation.gtf.gz gencode_exons.txt.gz` (GTF needs `gene_name` attributes) or omit `-e`. `deltapsi` = group 2 - group 1.

**LeafCutter2** (Buen Abad Najar 2025 *bioRxiv*) extends leafcutter with NMD-aware classification of unproductive splicing — useful when AS-NMD coupling is the question.
