# Which sort order a tool needs

No script; look the tool up, then sort with `routes/sort.md`.

| Operation | Required sort |
|-----------|---------------|
| `samtools index` | coordinate |
| `samtools fixmate -m` | name (`-n` or `-N`) or collate; coordinate input stops with `Coordinate sorted, require grouped/sorted by queryname` |
| `samtools markdup` | coordinate, after fixmate |
| Picard MarkDuplicates, ValidateSamFile | coordinate, or queryname in `-N` order |
| GATK MarkDuplicatesSpark | coordinate or queryname |
| `samtools mpileup`, `bcftools mpileup` | coordinate (name-sorted stops with `The input is not sorted`) |
| GATK HaplotypeCaller, Mutect2 | coordinate and indexed |
| htseq-count | `-r pos` for coordinate-sorted, `-r name` for name-sorted; the wrong `-r` gives different counts (5599 vs 2820 reads in one test) |
| featureCounts (unverified) | coordinate or name; `-p` for paired-end |
| umi_tools dedup | coordinate, indexed |
| fgbio GroupReadsByUmi | any order; template-coordinate avoids an internal re-sort (run `samtools fixmate -m` or `fgbio SetMateInformation` first) |
| fgbio CallMolecularConsensusReads | template-coordinate, as emitted by GroupReadsByUmi |
| Sniffles, cuteSV, Manta, Delly (unverified) | coordinate and indexed |
| Salmon alignment mode, RSEM (STAR `--quantMode TranscriptomeSAM`) (unverified) | mates adjacent, grouped by name; `-N` is the safe choice |

"Unverified" means not checked against the tool; confirm in its own documentation.

## collate vs sort -n

`collate` only makes mates adjacent (between-mate order undefined, `SO:unsorted GO:query`). Use it when the tool needs adjacency only: paired FASTQ extraction, re-aligning, fixmate, markdup pre-processing. It is not reliably faster than `sort -n` (0.74 s against 0.66 s on one slice); time both on your storage. Use `sort -n`/`-N` only when a real name order is required.

```bash
samtools collate -O -u in.bam tmp_prefix | \
    samtools fastq -1 R1.fq.gz -2 R2.fq.gz -0 /dev/null -s /dev/null -n -
```

Done when the tool's required order is named and the file is sorted that way.
