# Metrics and ENCODE 4 thresholds

Seven metrics, each with a route:

| Metric | Definition | Ideal | Acceptable | Reject | Route |
|--------|-----------|-------|------------|--------|-------|
| Nuclear reads (after dedup, no chrM) | Mapped reads (2 per fragment), MAPQ >= 30, non-chrM, deduped | >= 50M | 25-50M | < 25M | `routes/fragment-tss-frip.md` |
| Alignment rate | Mapped / total reads | >= 95% | 80-95% | < 80% | `routes/alignment-mt.md` |
| Mitochondrial fraction | chrM / total mapped | < 5% (Omni-ATAC), < 20% (standard) | 20-50% | > 50% | `routes/alignment-mt.md` |
| NRF | Distinct fragments / total fragments (single-end: 5' positions / reads) | >= 0.9 | 0.7-0.9 | < 0.7 | `routes/complexity.md` |
| PBC1 | Fragments seen once / distinct fragments | >= 0.9 | 0.7-0.9 | < 0.7 | `routes/complexity.md` |
| PBC2 | Fragments seen once / fragments seen twice | >= 3.0 | 1.0-3.0 | < 1.0 | `routes/complexity.md` |
| TSS enrichment (hg38) | Average signal at TSS / flanking | >= 7 | 5-7 | < 5 | `routes/tss-encode.md` |
| FRiP | Reads in peaks / total | >= 0.3 | 0.2-0.3 | < 0.2 | `routes/fragment-tss-frip.md` |
| Insert-size periodicity | NFR, mono-, di-nucleosome peaks visible | 3+ clear peaks | NFR + mono | Flat or single peak | `routes/fragment-tss-frip.md` |

Which bands are ENCODE and which are working convention:

- ENCODE 4 minimum depth is 25M single-end reads or 50M paired-end reads (25M fragments); the ideal/acceptable tiers here are this Skill's convention.
- ENCODE 4 sets no mt threshold; the bands above are working convention.
- NRF and PBC1 "ideal over 0.9" and PBC2 "ideal over 3" are ENCODE 4; the other bands are working convention.
- TSS enrichment cutoffs vary by annotation; ENCODE lists GRCh38 RefSeq. Mouse (mm10, GENCODE M21) at 5 or above is acceptable.
- ENCODE requires NFR and mononucleosome peaks in the insert-size distribution.

Methodology evolves: verify against the current ENCODE ATAC-seq standards before reporting. A full audit runs every route, then `routes/aggregate.md`.

Done when each metric has a value and a route has been run for it.
