---
name: bio-bam-statistics
description: Generate alignment statistics using samtools flagstat, stats, depth, coverage, and mosdepth. Use when assessing alignment quality, calculating coverage, or generating QC reports.
tool_type: cli
primary_tool: samtools
license: MIT
---

## Version Compatibility

Reference examples checked on: samtools 1.24, pysam 0.24.1, mosdepth 0.3.14, Picard 3.5.0, MultiQC 1.35 (written for pysam 0.22+, samtools 1.19+)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Install: `conda install -c bioconda samtools pysam mosdepth gnuplot perl-uri`. `plot-bamstats` ships with samtools and dies on `URI/Escape.pm` without `perl-uri`.

# BAM Statistics

**"Get alignment statistics and coverage from my BAM file"** -> Generate read counts, mapping rates, per-chromosome statistics, depth profiles, and coverage summaries.
- CLI: `samtools flagstat`, `samtools stats`, `samtools depth`, `samtools coverage` (samtools)
- Python: `pysam.AlignmentFile` with `pileup()` and `get_index_statistics()` (pysam)

Generate alignment statistics using samtools and pysam.

## Quick Summary Commands

| Question | Best tool | Why |
|----------|-----------|-----|
| Quick read counts by FLAG category | `samtools flagstat` | Fast; counts secondary+supp in totals |
| Per-chromosome counts | `samtools idxstats` | Fast (uses index); **counts secondary+supp** |
| Insert size, MAPQ, error, GC | `samtools stats` (`-r ref.fa` for GC-depth; CRAM needs `--reference ref.fa`) | Comprehensive; feeds MultiQC |
| Mean depth / breadth per contig or region | `samtools coverage` | One row per contig; uncovered bases count as zero |
| Per-position depth (small region) | `samtools depth -a` or pysam pileup | Slow on full genome |
| Per-position depth (genome-wide) | **`mosdepth`** | Multithreaded; preferred for whole genomes |
| Per-region coverage (BED) | `mosdepth --by regions.bed` | Production default |
| Coverage histogram / cumulative | `mosdepth -t 4 --no-per-base` | Single-pass histogram |
| Breadth at depth thresholds | `mosdepth --by regions.bed --thresholds 1,10,30,100` | Standard exome QC (`--thresholds` requires `--by`) |
| Targeted enrichment QC | `picard CollectHsMetrics` | PCT_OFF_BAIT, FOLD_80_BASE_PENALTY, AT/GC dropout |
| Cross-sample contamination / sample identity | `verifybamid2` / `somalier` | FREEMIX < 0.01 is the commonly used expectation |

### What Each Tool Counts (and Doesn't)

| Counting category | flagstat | stats | idxstats |
|-------------------|----------|-------|----------|
| Primary alignments | `primary` line (= `in total` minus secondary minus supp) | `raw total sequences` | not separable (mapped column includes secondary/supp) |
| Secondary | `secondary` line | filtered out | counted in mapped |
| Supplementary | `supplementary` line | filtered out | counted in mapped |
| Mapping rate denominator | `in total`, *including* secondary/supp | primary only | mapped+unmapped |

flagstat lines read `<QC-passed> + <QC-failed>`. Percentages use the QC-passed column; the QC-failed reads sit in the second column and are included in `in total`, `primary` and `raw total sequences`.

For long-read data where one read produces many supplementary alignments, the senior cross-check:
```
input_read_count = flagstat_total(passed + failed) - secondary - supplementary
                 = stats_raw_total_sequences
```
Reports of "the file has 1.2M reads" where the input was actually 800k with 400k supplementary chimeric splits trace to flagstat misinterpretation.

Mate-overlap handling differs between tools, so **mean depth differs 2x on the same paired-end BAM** (test BAM: 16.77x vs 8.86x). State which convention you report:

| Tool | Overlapping mates |
|------|-------------------|
| `samtools depth`, `samtools coverage`, pysam `pileup()`, `mosdepth --fast-mode` | counted twice (16.77x); no overlap option in `coverage` |
| `samtools depth -s`, `mosdepth` (default) | counted once (8.86x) |
| `samtools mpileup` | counted once at the default `-Q 13` (8.85x); `-x` disables it (16.75x). The lower-quality mate base is zeroed, so `-Q 0` defeats the removal (16.77x) |
| `bcftools mpileup` | counted once in `FORMAT/DP` (`-a FORMAT/DP`) and in the bases used for calling: 8.85x at the default, 16.77x with `-x` or `-Q 0`. `INFO/DP` is 16.77x either way (counted before overlap removal) |

## samtools flagstat

Fast summary of alignment flags.

```bash
samtools flagstat input.bam
```

Output:
```
10000000 + 0 in total (QC-passed reads + QC-failed reads)
9950000 + 0 primary
0 + 0 secondary
50000 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
9800000 + 0 mapped (98.00% : N/A)
9750000 + 0 primary mapped (97.99% : N/A)
9950000 + 0 paired in sequencing
4975000 + 0 read1
4975000 + 0 read2
9700000 + 0 properly paired (97.49% : N/A)
9720000 + 0 with itself and mate mapped
30000 + 0 singletons (0.30% : N/A)
15000 + 0 with mate mapped to a different chr
10000 + 0 with mate mapped to a different chr (mapQ>=5)
```
(samtools 1.13+ adds the `primary`, `primary duplicates`, and `primary mapped` lines shown above.) `properly paired` is the aligner-set proper-pair flag (both mates mapped, expected orientation and distance); `singletons` are mapped reads whose mate is unmapped; `duplicates` is non-zero only after a duplicate-marking step.

### Multi-threaded
```bash
samtools flagstat -@ 4 input.bam
```

### Machine-readable
`samtools flagstat -O tsv input.bam` prints `<passed>\t<failed>\t<label>` per line (percentages are separate `... %` rows).

### Summary table for many BAMs
```bash
printf 'Sample\tRecords\tQCfail\tPrimary\tPrimaryMapped\tProperPair\tPrimaryDup\n' > summary.tsv
for bam in *.bam; do
    sample=$(basename "$bam" .bam)
    samtools flagstat -O tsv "$bam" | awk -F'\t' -v s="$sample" '
        $3 ~ /^total/            {rec = $1 + $2; fail = $2}
        $3 == "primary"          {pri = $1 + $2}
        $3 == "primary mapped"   {pm = $1 + $2}
        $3 == "properly paired"  {pp = $1 + $2}
        $3 == "primary duplicates" {dup = $1 + $2}
        END {print s"\t"rec"\t"fail"\t"pri"\t"pm"\t"pp"\t"dup}' >> summary.tsv
done
```
`Records` is every alignment record (QC-passed + QC-failed, including secondary/supplementary); `Primary*` columns exclude secondary/supplementary and include QC-failed reads.

## samtools idxstats

Per-chromosome read counts. Uses the index; on 1.24 an unindexed BAM falls back to a slow full scan with a warning (pysam `get_index_statistics()` and mosdepth raise an error instead).

```bash
samtools idxstats input.bam
```

Output format: `chrom length mapped unmapped`
```
chr1    248956422    5000000    1000
chr2    242193529    4800000    800
chrM    16569        50000      100
*       0            0          150000
```

### Parse idxstats
Contig names differ between builds (`chrM` vs `MT`, `chrX` vs `X`), so check the names first; the recipes below exit 1 with a message instead of printing a false `0`.
```bash
# Total mapped alignments (includes secondary/supplementary)
samtools idxstats input.bam | awk '{sum += $3} END {print sum}'

# Mitochondrial percentage of mapped alignments (chrM or MT)
samtools idxstats input.bam | awk '
    $1 ~ /^(chr)?(M|MT)$/ {mt += $3; found = 1}
    {total += $3}
    END {if (!found) {print "no chrM/MT contig in idxstats" > "/dev/stderr"; exit 1}
         if (!total) {print "no mapped reads" > "/dev/stderr"; exit 1}
         printf "%.2f%% mitochondrial\n", mt/total*100}'

# Sex check (X/Y ratio; +1 avoids division by zero)
samtools idxstats input.bam | awk '
    $1 ~ /^(chr)?X$/ {x = $3; fx = 1}
    $1 ~ /^(chr)?Y$/ {y = $3; fy = 1}
    END {if (!fx || !fy) {print "no chrX/chrY contig in idxstats" > "/dev/stderr"; exit 1}
         printf "X:Y = %.2f\n", x/(y+1)}'
```

## samtools stats

Comprehensive statistics including insert size, base quality, and more.

```bash
samtools stats input.bam > stats.txt
```

### View Summary Numbers
```bash
samtools stats input.bam | grep "^SN" | cut -f 2-
```

Key summary fields:
- `raw total sequences` - Primary reads, QC-failed included (secondary/supplementary excluded)
- `reads mapped` - Mapped reads
- `reads mapped and paired` - Both mates mapped (NOT the proper-pair count)
- `reads properly paired` / `percentage of properly paired reads (%)` - Proper-pair flag
- `insert size average` / `insert size standard deviation` - Insert size mean and spread
- `inward oriented pairs` / `outward oriented pairs` / `pairs with other orientation` - Pair orientation
- `average length` - Mean read length
- `bases mapped (cigar)` - Mapped bases from the CIGAR (soft-clipped bases excluded; `bases mapped` includes them)
- `error rate` - Mismatch rate (mismatches / bases mapped (cigar))

### Generate Plots (with plot-bamstats)
```bash
samtools stats input.bam > stats.txt
plot-bamstats -p plots/ stats.txt
```

### Combine Samples (MultiQC)
```bash
multiqc . -o multiqc_out    # picks up samtools stats / flagstat / idxstats text outputs in the directory (checked on MultiQC 1.35)
```

### Stats for Specific Region
```bash
samtools stats input.bam chr1:1000000-2000000 > region_stats.txt
```

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
# When fragment length < 2 * read_length, R1 and R2 overlap (tool defaults: see the table above).
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
`mosdepth` excludes unmapped, secondary, QC-fail, and duplicate reads by default (`--flag 1796`); supplementary reads are NOT excluded (use `--flag 3844` to drop them too). Configurable via `--flag`. Memory ~ 4 bytes x longest chrom (1 GB for human chr1, 12+ GB for axolotl). Does not honor base quality; use `samtools depth -q INT` if needed. The summary `total` row covers only contigs that have reads (19000 of 22000 bp, 2.95x against 2.57x from `samtools coverage` / `depth -aa` on a BAM with one read-less contig); `--fast-mode` also counts deletion (D) bases as covered (69.97x vs 68.84x on an amplicon BAM). For a whole-reference mean use `samtools coverage` or `depth -aa`.

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

## pysam Python Alternative

Two inputs break a plain `AlignmentFile(path, 'rb')`: an unaligned BAM has no `@SQ` lines (`ValueError: file has no sequences defined`), so the read-only snippets pass `check_sq=False`; and a CRAM needs `reference_filename='ref.fa'`, without which iterating fails with `OSError: truncated file` (region and pileup calls also need a `.crai`; `samtools flagstat` needs no reference, `samtools stats` needs `--reference`). `examples/qc_report.py` takes the reference as an optional second argument and exits 1 with a message when a CRAM cannot be decoded.

### Count Reads

**Goal:** Reproduce `samtools flagstat` (QC-passed column) counts and rates.

**Approach:** Count secondary/supplementary records and QC-failed reads separately; all rates use QC-passed primary reads only, so they equal flagstat's `primary`, `primary mapped`, `properly paired` and `primary duplicates` lines. `examples/qc_report.py` is the full version.

```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb', check_sq=False) as bam:
    primary = mapped = paired = proper = 0
    for read in bam:
        if read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        primary += 1
        mapped += not read.is_unmapped
        paired += read.is_paired
        proper += read.is_proper_pair

if not primary:
    raise SystemExit('no QC-passed primary reads')
print(f'Primary QC-passed: {primary}')
print(f'Mapped: {mapped} ({mapped/primary*100:.2f}% of primary)')
if paired:
    print(f'Properly paired: {proper} ({proper/paired*100:.2f}% of primary paired)')
else:
    print('Single-end data: no paired reads')
```

### Per-Chromosome Counts
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    # needs an index (ValueError otherwise); counts include secondary/supplementary like idxstats
    for stat in bam.get_index_statistics():
        print(f'{stat.contig}: {stat.mapped} mapped, {stat.unmapped} unmapped')
```

### Depth in a Region (and at a Position)

**Goal:** Per-base depth, mean depth and breadth for a region that agree with `samtools depth -a`.

**Approach:** Fill a zero array over the region and record each pileup column; the mean and breadth divide by the region length, never by the number of covered columns.

**Reference (pysam 0.24.1, checked equal to `samtools depth -a` to 1e-9):**
```python
import pysam

def region_depth_stats(bam_path, chrom, start, end, thresholds=(10, 20), reference=None):
    """Depth over the 0-based half-open region [start, end); every base is in the denominator.
    A single position is the 1-bp region (pos - 1, pos). `reference` is the FASTA for a CRAM."""
    length = end - start
    if length <= 0:
        raise ValueError(f'empty region [{start}, {end})')
    depths = [0] * length
    with pysam.AlignmentFile(bam_path, 'rb', reference_filename=reference) as bam:
        # truncate=True: only columns inside the region (otherwise whole read footprints are returned)
        # max_depth: pysam's default cap is 8000; ignore_orphans=False: default drops paired reads
        #   that lack the proper-pair flag; min_base_quality=0: default is 13 (samtools depth: 0)
        # pileup.n counts deletions and ref-skips, so count real bases explicitly
        for col in bam.pileup(chrom, start, end, truncate=True, max_depth=1_000_000,
                              ignore_orphans=False, ignore_overlaps=False, min_base_quality=0):
            depths[col.reference_pos - start] = sum(1 for r in col.pileups if not r.is_del and not r.is_refskip)
    stats = {'length': length, 'covered': sum(d > 0 for d in depths),
             'mean_depth': sum(depths) / length, 'max_depth': max(depths)}
    stats['pct_covered'] = stats['covered'] / length * 100
    for t in thresholds:
        stats[f'pct_ge_{t}x'] = sum(d >= t for d in depths) / length * 100
    return stats

stats = region_depth_stats('input.bam', 'chr1', 1000000, 2000000)
print(f'Coverage: {stats["pct_covered"]:.1f}%  Mean depth: {stats["mean_depth"]:.1f}x  >=20x: {stats["pct_ge_20x"]:.1f}%')
```
Like `samtools depth`, this skips unmapped, secondary, QC-failed and duplicate reads and counts overlapping mates twice (the mean is 2x the `mosdepth` default on overlapping pairs). Python is slow beyond a few Mb: use `samtools coverage` or `mosdepth` for large regions.

### Insert Size Distribution

**Goal:** Compute the insert size distribution to assess library preparation quality.

**Approach:** Iterate QC-passed primary properly paired read1 records, accumulate template lengths into a Counter, then compute summary statistics.

**Reference (pysam 0.24.1):**
```python
import pysam
from collections import Counter

insert_sizes = Counter()

with pysam.AlignmentFile('input.bam', 'rb', check_sq=False) as bam:
    for read in bam:
        if read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if read.is_proper_pair and read.is_read1 and read.template_length > 0:
            insert_sizes[read.template_length] += 1

if not insert_sizes:
    raise SystemExit('no properly paired read1 records (single-end data, or proper-pair flag unset: see Insert Size Caveats)')
sizes = list(insert_sizes.keys())
mean_insert = sum(s * c for s, c in insert_sizes.items()) / sum(insert_sizes.values())
print(f'Mean insert size: {mean_insert:.0f}')
print(f'Min: {min(sizes)}, Max: {max(sizes)}')
```

## QC Thresholds Are Assay-Specific

A single "mapping rate > 95%" rule rejects valid ATAC, ChIP, RNA-seq, metagenomics, and aDNA samples. The threshold question is "is this rate normal for this assay?" not "is this rate above 95%?" Values below are literature ranges, not thresholds verified on this machine.

| Metric | WGS PCR-free | WGS PCR | WES | Targeted panel | Deep panel (UMI) | RNA-seq | scRNA (10x) | ATAC | ChIP | Long-read | aDNA |
|--------|--------------|---------|-----|----------------|------------------|---------|-------------|------|------|-----------|------|
| Mapping rate | >99% | >98% | >95% | >95% | >95% | >90% | >70% | >50% | >60% | >95% | 1-50% |
| Duplicate rate | <5% | 5-15% | 20-50% | 20-50% | 50-90% pre-consensus | (skip) | (use UMI) | 10-30% | 5-30% | n/a | 20-60% |
| Proper pair rate | >95% | >95% | >85% | >80% | >80% | >70% | n/a | >50% | >70% | n/a | >60% |
| Mean MAPQ | bimodal at 0/60 | bimodal | bimodal | bimodal | bimodal | bimodal incl 255 (STAR) | 0/1/3/255 | 30-55 | 30-55 | 30-50 | 20-40 |
| Mt fraction | 0.1-2% | 0.1-2% | <1% | <0.1% | <0.1% | varies | varies | **<10% (Omni-ATAC goal; original Buenrostro-2013 libraries were often majority-mito)** | <2% | n/a | varies |

Generic red flags whatever the assay: mapping rate < 80% (contamination, wrong reference, poor quality); `error rate` > 2%.

Mean MAPQ is misleading; the distribution is bimodal (0 and aligner-max). The fraction at MAPQ >= 30 is more informative:
```bash
samtools view -c -F 2308 -q 30 in.bam   # primary, mapped, MAPQ>=30
samtools view -c -F 2308 in.bam          # primary, mapped (denominator)
# For STAR/STARsolo, use -q 255 instead of -q 30 (255 is the unique-mapping sentinel)
```

## What Flagstat Does Not Reveal

A 99% flagstat mapping rate does NOT mean the data is usable. Common false-positive scenarios:

1. **Adapter readthrough**: short fragments (insert < 2 * read_length) sequence into adapter; aligners soft-clip the adapter portion and flag the read as MAPPED. `samtools stats` has no soft-clip field on 1.24 (a `grep "bases soft-clipped"` silently returns nothing), so count soft-clipped bases from the CIGAR of primary mapped reads (>5% is a rule of thumb for adapter contamination; local aligners and split reads also soft-clip):
   ```bash
   samtools view -F 2308 input.bam | awk -F'\t' '
       {n++; c = $6
        while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
            len = substr(c, 1, RLENGTH-1) + 0; op = substr(c, RLENGTH, 1); c = substr(c, RLENGTH+1)
            if (op ~ /[MIS=X]/) q += len; if (op == "S") s += len}}
       END {if (!n || !q) {print "no primary mapped reads: nothing to compute" > "/dev/stderr"; exit 1}
            printf "soft-clipped bases: %d of %d (%.2f%%)\n", s, q, s/q*100}'
   ```
   Denominator: query bases (M, I, S, =, X) of primary mapped records. Exits 1 with a message on an empty or unmapped BAM.
2. **Off-target enrichment** (capture/WES): `picard CollectHsMetrics` PCT_OFF_BAIT or PCT_SELECTED_BASES. Interval lists need the reference dictionary header (`picard BedToIntervalList I=targets.bed O=targets.interval_list SD=ref.dict`):
   ```bash
   picard CollectHsMetrics I=input.bam O=hs_metrics.txt R=ref.fa \
       BAIT_INTERVALS=baits.interval_list TARGET_INTERVALS=targets.interval_list
   ```
3. **Low-complexity pile-up**: telomere/centromere reads mass at MAPQ-0; counted as mapped but useless. Detect via MAPQ distribution.
4. **Cross-sample contamination**: `verifybamid2` estimates FREEMIX (> 1% degrades somatic calling; > 5% breaks germline calling are commonly cited cut-offs); `somalier` checks sample identity/relatedness (`extract` then `relate`, and has a `contamination` subcommand). Both need a whole-genome or exome BAM: on a small slice they report "No reads found in any of the regions" / too few markers.
   ```bash
   # --SVDPrefix must include the .dat suffix of the panel files (e.g. 1000g.phase3.10k.b38.vcf.gz.dat)
   verifybamid2 --SVDPrefix panel.vcf.gz.dat --Reference ref.fa --BamFile input.bam --Output sample_vb
   # FREEMIX is in sample_vb.selfSM

   somalier extract -s sites.vcf.gz -f ref.fa -d extracted/ input.bam    # FASTA must contain the sites' contigs
   somalier relate extracted/*.somalier
   ```
   (Flags checked against `--help` of verifybamid2 2.0.3 and somalier 0.3.5; not run end to end here because no whole-genome BAM was available.)
5. **Wrong reference build**: a BAM aligned to GRCh37 viewed against GRCh38 looks fine to flagstat but produces nonsense pileups. Compare `@SQ M5:` from BAM header with `samtools dict ref.fa` (if the header carries no M5, compare `@SQ` names and lengths with `ref.fa.fai`) -- see alignment-validation.

## Insert Size Caveats

`samtools stats` reports the IS section for every pair with both mates mapped and splits the pairs into `inward oriented`, `outward oriented` and `other orientation` counts (checked on a synthetic mate-pair library with the proper-pair flag set and unset: `insert size average` 2000.0, 100 outward pairs both times). So:
- Mate-pair libraries (RF orientation): IS is reported, with outward-oriented counts dominating. The pysam snippets and `qc_report.py` only look at properly paired reads, so they report nothing when the aligner leaves the proper-pair flag unset
- `qc_report.py` keeps insert sizes below `MAX_INSERT` = 8000 (the `samtools stats` default, `-i`); longer templates are dropped from its mean and median
- ATAC-seq: bimodal/multimodal expected (nucleosome ladder ~50/~180/~370 bp). Unimodal suggests poor transposition.
- RNA-seq: TLEN includes intron span -- mean meaningless
- Bisulfite (PBAT): orientation reversed; samtools may not flag proper pair

## Related Skills

- sam-bam-basics - View alignment files; aligner-aware MAPQ semantics
- alignment-indexing - idxstats requires index; secondary+supp counted
- alignment-validation - Insert size by library, contamination, sample-swap detection
- duplicate-handling - Library-aware duplicate rate expectations
- alignment-filtering - Filter before stats
- sequence-io/sequence-statistics - FASTA/FASTQ statistics
