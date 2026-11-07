# Sex check, spike-in, cell-cycle

Shown, not scripted.

Sex-chromosome check (sample swaps, clinical or biobank-scale work):

```bash
samtools idxstats sample.bam | awk '$1=="chrY"{print $3 / $2}'   # reads per bp
samtools view -c sample.bam chrX:73820651-73852753                # XIST locus, hg38
```

Female: chrY reads per bp near 0 and a high XIST count. Male: chrY at male coverage and a low XIST count. A discrepancy with the metadata flags a swap.

Spike-in (global accessibility shift expected, for example HDAC or DNMT inhibitors, differentiation): per-library normalization erases global shifts. Add Drosophila or E. coli chromatin before Tn5, align to a concatenated reference, count spike-in reads per sample and normalize by them. Workable spike-in fraction is 0.5-5% of reads; under 0.1% is unreliable and over 10% costs cellular reads.

Cell-cycle effect (proliferating lines such as K562, HEK293, HeLa): replication-associated opening shifts global accessibility by 5-15%. For bulk, score samples against a replication-origin atlas (Repli-seq peaks) and add the S-phase score to the DESeq2 design (`~Sphase + Condition`); for scATAC, regress it out as with Seurat CellCycleScoring.

Done when the answer states what the check showed and what it changes downstream.
