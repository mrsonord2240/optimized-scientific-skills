# NRF, PBC1, PBC2

```bash
python scripts/library_complexity.py sample.bam > complexity.json
```

Options: `--mode auto|paired|single`, `--min-mapq 30`, `--exclude-contigs chrM,MT,M,chrMT`. Input is the duplicate-retained BAM: flagged duplicates are counted, not skipped. Paired-end keys each proper pair by (chrom, fragment start, fragment end) from read 1; single-end keys by (chrom, 5' position, strand). The script prints JSON and exits 1 when no read passes the filters.

- Count fragments, not mate starts: per-mate counting undercounts complexity (0.35 by mate start against 0.62 by fragment on one unfiltered slice).
- PBC2 is `null` when no fragment is seen twice. That means the library is too sparse; PBC2 is unreliable below about 5M reads.
- NRF exactly 1.0 means the BAM was already deduplicated.
- To decide between re-prep and deeper sequencing, see `routes/preseq.md`.

Done when `complexity.json` has NRF, PBC1 and PBC2 from a pre-dedup BAM.
