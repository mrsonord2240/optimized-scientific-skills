## samtools depth

Per-position read depth.

### Basic Depth
```bash
samtools depth input.bam > depth.txt
```

Output: `chrom position depth`. Without `-a` only covered positions are printed, so **never average raw `samtools depth` output**: on the test BAM `awk '{s+=$3;n++} END{print s/n}'` gives 568x where the true mean is 16.8x.

### Depth at Specific Positions
```bash
samtools depth -r chr1:1000-2000 input.bam
```

### Include Zero-Depth Positions
```bash
samtools depth -a input.bam > depth_with_zeros.txt     # all positions of contigs that have reads
samtools depth -aa input.bam > depth_all_contigs.txt   # every position of every @SQ contig
```

### Mean Depth and Breadth
```bash
# Genome-wide mean and breadth: -aa puts every reference position in the denominator.
# For a region use `samtools depth -a -r chr1:1000-2000 input.bam` in the same pipe.
samtools depth -aa input.bam | awk '
    {s += $3; n++; if ($3 >= 10) c10++; if ($3 >= 20) c20++}
    END {if (!n) {print "no positions (empty BAM or no @SQ)" > "/dev/stderr"; exit 1}
         printf "mean %.2fx  >=10x %.2f%%  >=20x %.2f%%\n", s/n, c10/n*100, c20/n*100}'
```
For whole genomes this prints one line per base; use `samtools coverage` (`meandepth` column) or `mosdepth` (summary, `--thresholds`) instead.

### Maximum Depth Cap (Critical Trap)
```bash
# Default depth caps (checked on a 9500x stack; samtools 1.24, bcftools 1.24, pysam 0.24.1):
#   samtools depth         uncapped (-d/--max-depth is deprecated in 1.13+ and silently ignored)
#   samtools coverage      -d 1000000
#   mosdepth               uncapped
#   samtools mpileup       8000  -> samtools mpileup -d 1000000 -f ref.fa input.bam
#   bcftools mpileup       250   -> bcftools mpileup -d 1000000 -f ref.fa input.bam
#   pysam pileup()         max_depth=8000  -> pass max_depth=1_000_000
```
Silent capping under-reports deep positions (mean 850 instead of 1000 on the test stack).

Pipelines that break the mpileup/pysam cap: targeted oncology hotspots (5000-50000x), mitochondrial DNA (small genome, large read share), amplicon viral (ARTIC: 1000-100000x per amplicon), UMI-deduped capture (14000-17000x post-collapse), highly expressed transcripts (rRNA, mt-RNA).

### Overlapping Pair Correction
```bash
# When fragment length < 2 * read_length, R1 and R2 overlap (tool defaults: see the overlap table in SKILL.md "What Each Tool Counts").
# samtools depth counts the overlap twice; -s counts each template once:
samtools depth -s input.bam
```
Without `-s`, doubled support inflates somatic VAFs at sites covered by overlapping pairs (especially in fragmented samples: FFPE, cfDNA). `samtools mpileup` and `bcftools mpileup` remove overlaps by default; pass `-x` to disable (long form `--disable-overlap-removal` in samtools, `--ignore-overlaps` in bcftools).

### mosdepth (Modern Default)
```bash
mosdepth -t 4 sample input.bam                                                       # genome-wide per-base
mosdepth -t 4 --by exome.bed --thresholds 1,10,20,30,100 --no-per-base sample input.bam   # exome QC
mosdepth -t 4 --quantize 0:1:10:100: sample input.bam                                # CNV-style bands
mosdepth -t 4 -f ref.fa sample input.cram                                            # CRAM: needs the reference and a .crai (samtools index input.cram)
```
`mosdepth` excludes unmapped, secondary, QC-fail, and duplicate reads by default (`--flag 1796`); supplementary reads are NOT excluded (use `--flag 3844` to drop them too). Configurable via `--flag`. Memory ~ 4 bytes x longest chrom (1 GB for human chr1, 12+ GB for axolotl). Does not honor base quality; use `samtools depth -q INT` if needed. The summary `total` row covers only contigs that have reads (19000 of 22000 bp, 2.95x against 2.57x from `samtools coverage` / `depth -aa` on a BAM with one read-less contig); `--fast-mode` ignores internal CIGAR operations, so deletion (D) and spliced (N) bases count as covered (69.97x vs 68.84x on an amplicon BAM; 49.40x vs 17.67x on an RNA-seq BAM). For a whole-reference mean use `samtools coverage` or `depth -aa`.

### Depth from BED Regions
```bash
samtools depth -aa -b regions.bed input.bam   # every base of every region, zeros included
```
`-a` alone drops regions on contigs that have no reads (checked: 2810 of 3110 rows on a BAM with one read-less contig); `-aa` keeps them.

### Depth on Large Files
Restrict to a region, or subsample reproducibly (seed 42, 10%; approximate):
```bash
samtools depth -a -r chr1:1-10000000 input.bam                              # 10 Mb region
samtools view -s 42.1 -b -o sub.bam input.bam && samtools depth -a sub.bam
```

## samtools coverage

Per-chromosome or per-region coverage statistics.

```bash
samtools coverage input.bam
```

Output columns:
- `#rname` - Reference name
- `startpos` - Start position
- `endpos` - End position
- `numreads` - Number of reads
- `covbases` - Bases with coverage
- `coverage` - Percentage of bases covered
- `meandepth` - Mean depth (over the whole region, uncovered bases count as zero)
- `meanbaseq` - Mean base quality
- `meanmapq` - Mean mapping quality

### Coverage for Specific Region
```bash
samtools coverage -r chr1:1000000-2000000 input.bam
```

### Coverage from BED
`-b` is `--bam-list` in `samtools coverage` (not a BED option; a BED there fails with `Cannot open file list`). Loop over the BED (0-based start -> 1-based region) or use `mosdepth --by`:
```bash
while read -r chrom start end _; do
    case $chrom in ''|'#'*|track*|browser*) continue;; esac   # skip BED header lines
    samtools coverage -H -r "$chrom:$((start+1))-$end" input.bam
done < regions.bed
```

### Histogram Output
```bash
samtools coverage -m input.bam
```
