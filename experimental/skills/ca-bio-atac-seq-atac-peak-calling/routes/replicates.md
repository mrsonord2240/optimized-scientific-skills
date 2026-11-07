# Two-replicate IDR peaks

```bash
bash scripts/call_atac_peaks.sh rep1.bam rep2.bam EFFECTIVE_GENOME_SIZE blacklist.bed peaks_out chrom.sizes
```

The optional sixth argument creates a pooled bigWig; omit it if not needed. BAMs must be coordinate-sorted, indexed, deduplicated, MAPQ-filtered and mitochondrial-read-free. The script uses BAM read ends with `--shift -75 --extsize 150`, calls true and disjoint pseudo-replicates at p 0.01 and IDR 0.05, and writes blacklist-filtered conservative peaks and reproducibility ratios. It differs from the ENCODE tagAlign and summit pipeline; see `references/idr.md`. It does not apply the source's optional additional mean signal threshold or greylist. `MACS=macs2` and `IDR=/path/to/idr` override executables.

Done when `peaks_out/final/conservative.narrowPeak` and `peaks_out/idr/reproducibility.txt` exist and the ratios and filtering are reported.
