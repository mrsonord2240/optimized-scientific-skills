# Allele counts and Python pileups

```bash
python scripts/allele_counts.py input.bam chr1:1000000 > counts.txt
```

This indexed-BAM command accepts a 1-based position, filters MAPQ and base quality at 20, and does not apply BAQ. Functions in `scripts/pileup_helpers.py` take 0-based positions; see `references/pysam.md` for frequencies, region SNV candidates, text rendering and parameter equivalences. Reference skips and ambiguous bases are excluded from SNV counts.

Done when requested counts or position-level results are written with their filtering conventions.
