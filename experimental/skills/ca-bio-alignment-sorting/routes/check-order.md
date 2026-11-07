# Check sort order

```bash
python scripts/check_order.py input.bam [ref.fa]
```

Prints the `@HD` line, then walks every record and reports PASS or FAIL for coordinate order (reference order, then POS, unmapped last). It catches a mislabelled `SO:coordinate` BAM and reference blocks in the wrong order; the header alone does not. Exit 2 for a uBAM (no `@SQ`). For CRAM pass the reference or set `REF_PATH`, else pysam may report a misleading truncation error.

- Header only: `samtools view -H input.bam | grep '^@HD'`. `SO:unsorted` also covers tag and template-coordinate sorts.
- Quick indexability test: `samtools index input.bam /tmp/check.bai && echo indexable` fails with `Unsorted positions on sequence ...` within a contig, but accepts contigs in a different order from `@SQ`; use the script for the full check.
- Truncated or empty output after a long sort: `samtools quickcheck -v sorted.bam` (rc 16 truncated, rc 4 zero-byte).
- A name-sorted file's order is in `SS:` of the `@HD` line (`natural` or `lexicographical`).

Done when the script printed PASS, or you re-sorted with `routes/sort.md` after a FAIL.
