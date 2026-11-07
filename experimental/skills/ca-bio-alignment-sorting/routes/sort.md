# Sort a BAM or CRAM

```bash
bash scripts/sort_bam.sh coordinate in.bam sorted.bam [threads] [mem_per_thread]
```

Orders: `coordinate` (indexing, visualization, variant calling), `name-natural` (`-n`, samtools-only consumers), `name-ascii` (`-N`, Picard, GATK, htsjdk), `tag:CB` (group by a tag value, e.g. 10x barcode), `template-coordinate` (fgbio GroupReadsByUmi). Threads default to 4 and `-m` to 2G. The script sorts to a temp file, runs `samtools quickcheck -v`, publishes, prints the `@HD` line, and indexes a coordinate sort.

| Order | Header written |
|-------|----------------|
| coordinate | `SO:coordinate` |
| name-natural | `SO:queryname SS:queryname:natural` |
| name-ascii | `SO:queryname SS:queryname:lexicographical` |
| tag | `SO:unsorted SS:unsorted:TAG:coordinate` |
| template-coordinate | `SO:unsorted SS:unsorted:template-coordinate GO:query` |

- Picard 3.5.0 `MarkDuplicates` on `sort -n` output fails with `Alignments added out of order ... Sort order is queryname` whenever natural and ASCII orders differ; `-N` passes.
- `-@ N` adds N worker threads; output records are identical for any `-@`. `-m` is a per-thread budget; data beyond it spills to temp files and the result is unchanged. Values below `1M` are rejected.
- `SORT_TMP=<prefix>` puts temp files on fast local storage. `-T` is a prefix; an existing directory gets its temp files inside it.
- CRAM output: name the output `*.cram` and set `REF=ref.fa`. The script does not index a CRAM.
- `template-coordinate` runs `collate | fixmate -m` first because the sort needs `MC` tags (`no MC tag. Please run samtools fixmate on file first.` otherwise).
- Sort and index in one pass with a bare `samtools sort --write-index -o sorted.bam in.bam`; it writes `.csi`, not `.bai`.

Done when `quickcheck` passed, the printed `@HD` line shows the order you asked for, and a coordinate sort has a `.bai`.
