## leafcutter Junction Quantification

**Goal:** Detect junctions and intron clusters annotation-free for downstream cluster-level usage.

**Approach:** Extract junctions per BAM with regtools, write filenames into a list, then cluster introns sharing splice sites.

**Prerequisite:** the BAMs must carry an `XS` strand tag (STAR: `--outSAMstrandField intronMotif`; HISAT2 writes XS by default). Without it `regtools -s XS` records every junction with strand `?`, leafcutter drops them all, and the run still exits 0 with zero clusters. Check before clustering: `tail -n +2 s.junc | cut -f6 | sort | uniq -c` must show `+` and `-`, not only `?`. `leafcutter_cluster_regtools.py` is in the `clustering/` directory of the leafcutter repository (github.com/davidaknowles/leafcutter).

```bash
for bam in *.bam; do
    regtools junctions extract -a 8 -m 50 -s XS "$bam" -o "${bam%.bam}.junc"
done
ls *.junc > juncfiles.txt

python leafcutter_cluster_regtools.py \
    -j juncfiles.txt \
    -o leafcutter \
    -m 50 \
    -l 500000
```

`-a 8` = 8nt anchor minimum (raise to 12 for stricter; lower to 6 for microexon-friendly). `-m 50` = minimum junction reads per cluster (on small test data it can leave zero clusters; `-m 5` gave 115 on four chrX BAMs). `-l 500000` = max intron length (relevant for long brain-gene introns; raise for genes like DSCAM, ROBO2, ANK3).
