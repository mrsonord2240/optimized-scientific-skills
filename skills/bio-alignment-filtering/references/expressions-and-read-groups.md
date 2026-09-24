## Expression Filtering

`samtools view -e EXPR` (or `--expr`, since samtools 1.12) supports arbitrary expression filtering on tags, FLAG, MAPQ, RNAME, CIGAR, etc. Powerful for filtering by `NM`, `AS`, `NH`, `cs`, etc. that the FLAG-based filters cannot reach (the `sclen` keyword used below is documented from samtools 1.16):
```bash
# Reads with >=2 mismatches (NM tag)
samtools view -e '[NM] >= 2' in.bam

# Soft clip on the left, on chr1
samtools view -e 'cigar=~"^[0-9]+S" && rname=="chr1"' in.bam

# Combine with FLAG and MAPQ
samtools view -F 2308 -q 30 -e '[NM] <= 5 && [AS] >= 100' in.bam

# Drop reads with low mapped fraction (samtools-internal helpers)
samtools view -e 'sclen / qlen < 0.2' in.bam
```

Note: `![NM]` is true only if NM is missing (checked on 1.24); NULL values from missing tags propagate through arithmetic. A tag the file never carries (`[XY]`, or `[NH]` on BWA output) matches nothing and gives an empty result with exit code 0 and no warning.

## Filter by Read Group
`-r` takes the read group **ID** (the `ID:` of an `@RG` header line), not a library name, and also outputs reads that carry no RG tag:
```bash
samtools view -r SRR702039 in.bam              # single read group ID (plus untagged reads)
samtools view -e '[RG]=="SRR702039"' in.bam    # strictly that read group
samtools view -R rg_list.txt in.bam            # multiple IDs via file (one ID per line)
samtools view -l LIBRARY in.bam                # by library (@RG LB: field)
```
Samtools 1.23 adds `-n` (`--exclude-no-read-group`) to drop untagged reads when `-r`/`-R` is used (the 1.24 `--help` spells the long option `--exclude-no-read_group`; the hyphenated form also works on 1.24).
