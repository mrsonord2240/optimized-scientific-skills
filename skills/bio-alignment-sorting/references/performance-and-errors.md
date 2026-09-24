# Performance and Common Errors

Read this reference when choosing samtools compression or diagnosing a sort failure.

## Compression Level Decision

| Level | Use | Wall-time vs default | Size vs default |
|-------|-----|----------------------|------------------|
| `-l 0` / `-u` | Pipe between samtools tools | often fastest on local storage | many times larger (39x on the local test slice) |
| `-l 1` | Final output if disk is cheap | faster | ~8% larger |
| `-l 6` | Default | baseline | baseline |
| `-l 9` | Archival, write-once | ~10x slower | ~4% smaller |

Measured once on a 192k-read 1000 Genomes slice on local WSL storage with samtools 1.24 (`-@ 0`); data and storage change the ratios, so treat them as direction, not promises.

```bash
# WRONG -- pipe re-compresses then decompresses every step
samtools fixmate -m in.bam - | samtools sort -o out.bam

# RIGHT -- uncompressed (-u) between piped samtools commands
samtools fixmate -m -u in.bam - | samtools sort -o out.bam
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `out of memory` | Insufficient RAM | Use `-m` to limit per-thread memory; measure peak memory on the target host |
| `-m setting ... is less than the minimum required (1M)` | `-m` too small | Use `-m 1M` or more |
| `disk full` | Temp files filling disk | Use `-T` to point the temp prefix at a different disk |
| `truncated file` | Interrupted sort | Re-run sort from the original (keep it until the output passes `samtools quickcheck -v`) |
| `Unsorted positions on sequence` (index) | Not coordinate-sorted, or `@HD` says so wrongly | Re-sort; verify records, not the header |
| `NO_COOR reads not in a single block` / `cannot be indexed` | Indexing a `-n`, `-N` or `-t` output | Sort by coordinate first |
| `Alignments added out of order` (Picard) | Name-sorted with `-n` (natural) | Re-sort with `-N` |
| `no MC tag. Please run samtools fixmate` | `--template-coordinate` without mate tags | `samtools fixmate -m` first |
