# Fetch reference sequences

```bash
samtools faidx reference.fa chr1:1000-2000 > region.fa
```

Coordinates are 1-based inclusive; pysam `FastaFile.fetch` is 0-based half-open. Use `-i` for reverse complement, or `-i --mark-strand no` to keep the header plain. Pass multiple regions or contig names to extract a subset. Index the subset and create its dictionary if required. Verify coordinates against the FAI: an entirely out-of-range request can return an empty record with exit 0. Use `cut -f1,2 reference.fa.fai > chrom.sizes` for sizes.

Done when requested sequences are nonempty and match the requested coordinates.
