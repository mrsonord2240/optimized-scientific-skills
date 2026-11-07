# Errors, speed and size

| Error | Cause | Fix |
|-------|-------|-----|
| `out of memory` | `-m` times threads too high | Lower `-m`; measure peak memory on your host |
| `-m setting ... is less than the minimum required (1M)` | `-m` too small | `-m 1M` or more |
| `disk full` | Temp files filling the disk | `SORT_TMP=` on another disk |
| `truncated file` | Interrupted sort | Re-run from the original; check with `samtools quickcheck -v` |
| `Unsorted positions on sequence` (index) | Not coordinate-sorted, or `@HD` wrong | Re-sort; verify records (`routes/check-order.md`) |
| `NO_COOR reads not in a single block` / cannot be indexed | Indexing a `-n`, `-N` or `-t` output | Sort by coordinate first |
| `Alignments added out of order` (Picard) | Name-sorted with `-n` | Re-sort with `-N` |
| `no MC tag. Please run samtools fixmate` | `--template-coordinate` without mate tags | `samtools fixmate -m` first (the script does) |

Slow sort: raise `-@` and `-m`, put temp files on fast local disk.

Compression of final output (one 192k-read slice, local storage; direction, not promise):

| Level | Use | Time vs default | Size vs default |
|-------|-----|-----------------|-----------------|
| `-l 0` / `-u` | Pipe between samtools tools only | often fastest | 39x larger |
| `-l 1` | Final output if disk is cheap | faster | about 8% larger |
| `-l 6` | Default | baseline | baseline |
| `-l 9` | Archival | about 10x slower | about 4% smaller |

```bash
# WRONG: recompresses at every step
samtools fixmate -m in.bam - | samtools sort -o out.bam
# RIGHT
samtools fixmate -m -u in.bam - | samtools sort -o out.bam
```

Done when the failing command has run to completion.
