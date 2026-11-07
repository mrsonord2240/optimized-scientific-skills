# Paired 16S reads to an ASV table

```bash
Rscript scripts/dada2_asv.R --trunc-len 220,200 --out dada2_out run1_trimmed run2_trimmed
```

One directory per sequencing run, each holding that run's primer-trimmed `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`. `--trunc-len` is required; set it with `routes/truncation.md` (220,200 suits V4 on 2x250 after primer removal). Per run the script does `filterAndTrim` (`maxEE` 2,2, `truncQ` 2, `maxN` 0, PhiX removed), `learnErrors` for forward and reverse separately, `dada`, `mergePairs`; it then `mergeSequenceTables` across runs and runs one `removeBimeraDenovo` (consensus).

Options: `--max-ee 2,5` (loosen the reverse for long amplicons), `--amplicon-len N` (refuses a truncLen sum below N+12), `--pool pseudo` (cross-sample singletons; default FALSE misses them), `--pattern-f/--pattern-r` for other file suffixes, `--min-len`. It also stops if truncLen exceeds the longest primer-trimmed read.

Writes `seqtab_nochim.rds` (samples x ASVs), `read_tracking.csv`, `error_fit_*.png`, `asv_lengths.txt`.

- Open `error_fit_*.png`: observed points must track the fitted line and fall with quality. On NovaSeq, NextSeq or iSeq (about 4 quality bins) the loess fit can go non-monotonic; enforce a monotonic error matrix (nf-core/ampliseq `--illumina_novaseq`, or set sub-max-Q entries to the max-Q error) before trusting the denoising.
- Read `read_tracking.csv`: a near-zero `merged` column is a truncLen budget problem, not low diversity.
- Chimeras are many ASVs but few reads; retained reads normally 0.8-0.99. Much lower means primers are still on. The script prints the fraction.
- `removeBimeraDenovo` is not exhaustive: in a planted-truth check (DADA2 1.34.0) a low-abundance chimera (0.6% of reads) survived consensus, pooled and per-sample alike. Do not call the table chimera-free; inspect low-abundance ASVs that splice two abundant ones.
- Carry run forward as a batch covariate.
- Copy-number correction: report, do not assume.
- Install: `BiocManager::install(c('dada2', 'decontam'))`, plus `ggplot2`.

Done when `seqtab_nochim.rds` exists, the error plots were inspected, and `merged` in `read_tracking.csv` is sensible.
