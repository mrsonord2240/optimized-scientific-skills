# Errors

Messages are verbatim from samtools 1.24. Each stops the tool with exit 1 and leaves a partial output file; delete it before re-running.

| Error | Cause | Fix |
|-------|-------|-----|
| `[bam_mating_core] ERROR: Coordinate sorted, require grouped/sorted by queryname` | fixmate input is coordinate-sorted | `samtools sort -n` or `collate` first |
| `samtools markdup: error, no ms score tag. Please run samtools fixmate on file first.` | fixmate ran without `-m`, or a re-sort dropped the tags | Re-run `samtools fixmate -m` |
| `samtools markdup: error, no MC tag. Please run samtools fixmate on file first.` | Mate CIGAR tag lost | Re-run `samtools fixmate -m` |
| `samtools markdup: error, queryname sorted, must be sorted by coordinate.` | markdup input still name-sorted | `samtools sort` after fixmate |
| `Cannot open intermediate file "tmpdir/collate.0000.bam"` | Temp directory for `collate` or `sort -T` does not exist | `mkdir -p tmpdir` |
| `Mate mapping quality (MQ) tag not present` (fgbio) | GroupReadsByUmi input has no MQ tag | `samtools fixmate -m` or `fgbio SetMateInformation` first |
| `fetch called on bamfile without index` (umi_tools) | Input not sorted and indexed | `samtools sort`, then `samtools index` |
| `NullPointerException ... getReadGroupId()` (Picard) | No `@RG` line | `samtools addreplacerg` |

Done when the command that failed has run to completion.
