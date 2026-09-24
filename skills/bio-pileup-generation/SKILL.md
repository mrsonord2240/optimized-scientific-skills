---
name: bio-pileup-generation
category: Data Analysis
description: Generate pileup data for variant calling using samtools mpileup and pysam. Use when preparing data for variant calling, analyzing per-position read data, or calculating allele frequencies.
tool_type: cli
primary_tool: samtools
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples checked on: samtools 1.24, bcftools 1.24, pysam 0.24.1 (earlier tested with 1.19+ / 0.22+)

Install: `conda install -c bioconda samtools bcftools`, `pip install pysam`. The reference must be indexed (`samtools faidx reference.fa`) and the BAM must be indexed for `-r` and for `bam.pileup()`.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Pileup Generation

Generate pileup data for variant calling and position-level analysis.

**"Generate pileup from BAM"** -> Produce per-position read summaries showing depth, bases, and qualities.
- CLI: `samtools mpileup -f ref.fa input.bam`
- Python: `bam.pileup(chrom, start, end)` (pysam)

**"Count alleles at a position"** -> Extract per-base read support at a specific genomic coordinate.
- Python: iterate `pileup_column.pileups` and count bases (pysam)

## What is Pileup?

Pileup shows all reads covering each position in the reference, used for:
- Variant calling (with bcftools)
- Coverage analysis
- Allele frequency calculation
- SNP/indel detection

## samtools mpileup vs bcftools mpileup (Deprecation)

`samtools mpileup -g/-u` (BCF output for variant calling) was **deprecated in samtools 1.9 and removed in 1.15** (the option no longer exists; the usage/manpage directs users to `bcftools mpileup`) -- the genotype-likelihood code now lives in `bcftools mpileup`, which keeps mpileup logic versioned alongside `bcftools call` and avoids version-skew bugs.

| Use case | Recommended tool |
|----------|------------------|
| Quick allele counts at known sites | `samtools mpileup` or pysam pileup |
| Germline variant calling (small genomes, simple cohorts) | `bcftools mpileup` -> `bcftools call` |
| Germline WGS / WES production | DeepVariant or HaplotypeCaller (not mpileup) |
| Somatic SNV/indel | Mutect2 / VarDict / VarScan2 (direct from BAM) |
| Long-read small variants | clair3 / DeepVariant ONT (direct from BAM) |
| Long-read SV | Sniffles / cuteSV (direct from BAM) |
| Ultra-low-frequency (ctDNA / MRD) | fgbio consensus -> `bcftools call` or hot-spot Mutect2 |
| Per-position allele counts (custom) | pysam pileup |

`samtools mpileup` (without `-g`) is still the standard tool for human-readable per-position read summaries. Its output is **text, not VCF/BCF**, so it cannot feed `bcftools call` (see Maximum Depth below).

### Basic Pileup
```bash
samtools mpileup -f reference.fa input.bam > pileup.txt
```

### Pileup Specific Region
```bash
samtools mpileup -f reference.fa -r chr1:1000000-2000000 input.bam
```

### Regions from BED
```bash
samtools mpileup -f reference.fa -l targets.bed input.bam
```

### Multiple BAM Files
```bash
samtools mpileup -f reference.fa sample1.bam sample2.bam sample3.bam > pileup.txt
```

## Output Format

Text pileup format: 3 shared columns, then 3 columns per input BAM (6 columns for one BAM, 9 for two):
```
chr20   1400300  G    6    ,,,.,,    BBB<BB
chr20   1400301  C    6    ,,,.,,    BB<BBB
```

| Column | Description |
|--------|-------------|
| 1 | Chromosome |
| 2 | Position (1-based) |
| 3 | Reference base (`N` if `-f` is missing or lacks the contig; lower case where the FASTA is soft-masked) |
| 4 | Read depth (= number of read symbols in column 5 = number of characters in column 6) |
| 5 | Read bases |
| 6 | Base qualities (ASCII - 33; BAQ-adjusted unless `-B`) |

A covered position whose bases are all filtered prints depth 0 and `*` in columns 5 and 6 (e.g. `chr22 1952 T 0 * *`); uncovered positions are absent unless `-a`/`-aa`.

### Read Bases Encoding

| Symbol | Meaning |
|--------|---------|
| `.` | Match on forward strand |
| `,` | Match on reverse strand |
| `ACGT` | Mismatch (uppercase = forward) |
| `acgt` | Mismatch (lowercase = reverse) |
| `^Q` | Start of read (Q = MAPQ + 33 as ASCII, e.g. `^]` = MAPQ 60) |
| `$` | End of read (follows the read's last base) |
| `+2AC` | Insertion of the 2 bases AC after this base (adds no depth; lowercase = reverse) |
| `-3CGT` | Deletion of the 3 reference bases CGT after this base |
| `*` | Deleted base (counted in depth; `#` with `--reverse-del` on the reverse strand) |
| `>` / `<` | Reference skip (intron, `N` in CIGAR), forward / reverse read; counted in depth |

## Quality Filtering Options

The defaults also drop data silently; they are listed in Pileup Options and Defaults below.

### Minimum Mapping Quality
```bash
samtools mpileup -f reference.fa -q 20 input.bam
```

### Minimum Base Quality
```bash
samtools mpileup -f reference.fa -Q 20 input.bam
```

### Combined Quality Filters
```bash
samtools mpileup -f reference.fa -q 20 -Q 20 input.bam
```

### Maximum Depth (Critical Trap)
```bash
# The default -d of each tool (Pileup Options and Defaults) silently truncates targeted / mt-DNA / amplicon /
# UMI-deduped data; bcftools' is far lower. Set -d explicitly in whichever you use
samtools mpileup -f reference.fa -d 0 input.bam        # no cap
samtools mpileup -f reference.fa -d 1000000 input.bam  # explicit high cap

# WRONG -- samtools mpileup writes text; bcftools call reads VCF/BCF and stops with
# "Failed to read from standard input: unknown file type"
samtools mpileup -f ref.fa in.bam | bcftools call -mv

# RIGHT -- bcftools mpileup writes the BCF that bcftools call reads; set -d explicitly
bcftools mpileup -d 1000000 -f ref.fa in.bam | bcftools call -mv
```
pysam `pileup()` also caps at `max_depth=8000`, and `max_depth=0` is **not** unlimited (still 8000); pass a large number such as `max_depth=1000000`.

## BAQ: Base Alignment Quality (Critical Default)

When `-f ref.fa` is passed, BAQ is enabled by default. BAQ Phred-scales the probability that a base is misaligned (HMM realignment against the reference over a small window) and reduces base quality near indels. Tradeoffs: slower; suppresses FP SNVs near indels; hurts indel detection sensitivity. In a **text** pileup it can hide an indel completely: with BAQ on, the bases of reads carrying a 3 bp deletion fall below `-Q 13` next to it (depth 8 -> 4, no `-3CGT` marker); `-B` shows all 8 reads and the marker.

| Flag | Behavior |
|------|----------|
| (default with `-f`) | BAQ on (computed against the reference; an existing BQ tag is reused) |
| `-B` / `--no-BAQ` | Disable BAQ -- raw qualities |
| `-E` / `--redo-BAQ` | Force recompute, ignoring existing BQ tags (cannot be combined with `-B`) |

**BAQ ON for:** short-read germline SNV (BWA, Bowtie2, HISAT2), short-read somatic SNV.

**BAQ OFF (`-B`) for:** long-read variant calling (ONT, PacBio HiFi), SV calling, RNA-seq near splice junctions, viral / amplicon, ultra-deep ctDNA from consensus reads (consensus quality already inflated), aDNA (qualities pre-rescaled by mapDamage).

### Library-Typed Flags Cheat Sheet

| Library | Flags |
|---------|-------|
| Short-read germline WGS (BWA) | `-q 20 -Q 20 -d 0` (BAQ on default) |
| Short-read tumor WGS | `-q 1 -Q 13 -d 0 -B` (low MAPQ kept; BAQ off) |
| Amplicon viral (ARTIC) | `-aa -A -d 600000 -B -Q 20` |
| Capture / exome | `-q 20 -Q 20 -d 0` |
| Long-read ONT R10.4+ | `-q 30 -Q 0 -B -d 0`; for `bcftools mpileup` add `--max-BQ 30` (its `ont` preset value) |
| PacBio HiFi | `-q 20 -Q 0 -B -d 0` |
| RNA-seq variants | `-q 20 -Q 20 -B -d 0` |
| Forensic / aDNA | `-q 0 -Q 0 -A -d 0 -B` |

## Variant Calling Pipeline (Modern: bcftools mpileup)

**Goal:** Call variants from alignment data using the pileup-based approach.

**Approach:** Use `bcftools mpileup` (not `samtools mpileup -g`) so genotype-likelihood code is co-versioned with `bcftools call`. Apply quality and depth caps explicitly; annotate FORMAT fields needed for downstream filtering.

### Modern Germline Calling
```bash
bcftools mpileup -f reference.fa -d 1000000 -q 20 -Q 20 \
    --annotate FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
    input.bam | \
  bcftools call -mv -Oz -o variants.vcf.gz
bcftools index -t variants.vcf.gz
```

### BCF Intermediate (re-run `bcftools call` without repeating mpileup)
```bash
bcftools mpileup -f reference.fa -d 1000000 -Ob -o raw.bcf input.bam
bcftools call -mv raw.bcf -o variants.vcf
```

### Multi-Sample Joint Calling
```bash
bcftools mpileup -f reference.fa --threads 4 -d 100000 -q 20 -Q 20 \
    -a FORMAT/AD,FORMAT/DP s1.bam s2.bam s3.bam | \
  bcftools call -mv --threads 4 -Oz -o joint.vcf.gz
```

Joint calling over all samples beats merging per-sample VCFs. `-d` is per file: set it above the expected depth, but `-d 1000000` with several samples makes bcftools warn "Potential memory hog" (checked with 2 samples).

### Parallel by Contig
```bash
# Merge in header order: a chr*.vcf.gz glob sorts chr1, chr10, chr11, chr2 ... and bcftools concat
# accepts that silently (exit 0, valid index). Keep contigs.txt to primary contigs (no ':' or '*' names).
samtools idxstats input.bam | cut -f1 | grep -v '^\*$' > contigs.txt
xargs -a contigs.txt -P 4 -I{} bash -c \
    'set -o pipefail; bcftools mpileup -f reference.fa -r {} -d 1000000 input.bam | bcftools call -mv -Oz -o {}.vcf.gz' \
  && sed 's/$/.vcf.gz/' contigs.txt > vcf_list.txt \
  && bcftools concat -f vcf_list.txt -Oz -o all.vcf.gz && bcftools index all.vcf.gz
```

For somatic / low-VAF, prefer Mutect2 / Strelka2 / DeepVariant -- materially better than mpileup-based callers.

### Overlap Detection Defaults

When fragment length < 2 * read_length, R1 and R2 overlap. Both `samtools mpileup` and `bcftools mpileup` enable overlap detection by default and count overlapping bases once; pass `-x` to disable (long forms: `--disable-overlap-removal` / `--ignore-overlaps-removal` in samtools 1.24, `--ignore-overlaps` in bcftools). Disabling overlap correction can inflate somatic VAFs at sites covered by overlapping pairs (especially in cfDNA / FFPE).

## pysam Python Alternative

`bam.pileup(contig, start, stop)` is 0-based, half-open; pass `truncate=True` or it also yields the columns of every overlapping read outside the region. `pileup_column.pos` is 0-based.

**Depth is `len(pileup_column.pileups)`, not `pileup_column.n`.** `n` counts reads before the base-quality filter and overlap removal, so it exceeds the mpileup depth (on the human test BAM it differs at 1087 of 1157 positions). `len(pileup_column.pileups)` (= `get_num_aligned()`) equals the `samtools mpileup` depth column when the parameters below match.

**BAQ is switched on by `fastafile`, and only under pysam's default stepper.** The default stepper is `'samtools'` (not `'all'`); leave `stepper` unset or pass `'samtools'`. Without `fastafile`, pysam defaults equal `samtools mpileup -B`. With `fastafile=pysam.FastaFile(ref)` BAQ is applied and the pileup equals `samtools mpileup -f ref.fa`; add `compute_baq=False` to get `-B` back. An explicit `stepper='all'` applies no BAQ even with `fastafile` (it equals `-B`) and also drops overlap and orphan handling, so it matches neither mpileup default. Checked position by position (depth and base counts) on 6 BAMs (human DNA, spliced RNA-seq, 1000G, ARTIC nanopore, 2 synthetic): 0 differing positions for each mapping in the table.

| samtools mpileup | `bam.pileup()` argument | Note |
|------------------|-------------------------|------|
| `-f ref.fa` (BAQ on) | `fastafile=FastaFile(ref)` | default stepper (`'samtools'`) only; `stepper='all'` gives no BAQ |
| `-B` | no `fastafile`, or `compute_baq=False` | |
| `-E` | `fastafile=FastaFile(ref), redo_baq=True` | recomputes over an existing BQ tag; without it the tag is reused |
| `-Q 13` (default) | `min_base_quality=13` (default) | pysam default matches |
| `-q N` | `min_mapping_quality=N` | |
| `-d N` (default 8000) | `max_depth=N` (default 8000) | `0` is not unlimited (see Maximum Depth) |
| `-x` | `ignore_overlaps=False` | default `True` matches |
| `-A` | `ignore_orphans=False` | default `True` matches |
| `--ff` | `flag_filter=INT` (default 1796 = the same four flags) | `--ff 0` = `flag_filter=0`; `stepper='nofilter'` also drops the orphan filter, so it is not `--ff 0` |
| `--rf` | `flag_require=INT` | any of the bits set, like `--rf` |
| `-C 50` | `adjust_capq_threshold=50` | |

Each read in `pileup_column.pileups` needs three checks in this order: **`is_refskip` first**, because pysam sets `is_del=True` on spliced-read `N` skips as well; then `is_del`; else the base at `query_position`. Testing `is_del` first reports intron positions as deletions. `pileup_read.indel` is the length of the insertion (>0) or deletion (<0) that follows the base.

### Basic Pileup
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1001000, truncate=True):
        print(f'{pileup_column.reference_name}:{pileup_column.pos + 1} depth={len(pileup_column.pileups)}')
```

### Access Reads at Position
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1000001, truncate=True):
        print(f'Position: {pileup_column.pos}')
        print(f'Depth: {len(pileup_column.pileups)}')

        for pileup_read in pileup_column.pileups:
            if pileup_read.is_refskip:
                print('  Reference skip')
            elif pileup_read.is_del:
                print('  Deletion')
            else:
                aln = pileup_read.alignment
                qpos = pileup_read.query_position
                strand = '-' if aln.is_reverse else '+'
                print(f'  {aln.query_name} {strand} {aln.query_sequence[qpos]} (Q{aln.query_qualities[qpos]})')
```

### Shipped pysam helpers

`examples/pileup_helpers.py` contains the reusable `allele_counts`,
`allele_frequency`, `find_variants`, and `pileup_text` functions. For
`allele_counts`, use a **0-based pos: pass 1-based position minus 1**;
`find_variants` returns 1-based VCF-style positions. Read-base
`N`, reference skips, and reference-`N` sites are excluded from the
SNV-oriented counts, so frequency denominators agree across the helpers.

```python
from pileup_helpers import allele_counts, allele_frequency, find_variants

counts = allele_counts('input.bam', 'chr1', 1000000 - 1, min_base_quality=20)
frequency = allele_frequency('input.bam', 'chr1', 1000000 - 1, min_base_quality=20)
calls = find_variants('input.bam', 'reference.fa', 'chr1', 999000, 1000000)
```

Run the command-line counter from its directory with
`python allele_counts.py input.bam chr1:1000000` (1-based, MAPQ/baseQ >= 20).
Run `python self_test.py` there to build a temporary BAM and verify all helpers.

### Pileup with Quality Filtering
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1001000,
                                     truncate=True,
                                     min_mapping_quality=20,
                                     min_base_quality=20):
        print(f'{pileup_column.pos + 1}: {len(pileup_column.pileups)}')
```

### Generate Pileup Text

```python
from pileup_helpers import pileup_text

for row in pileup_text('input.bam', 'reference.fa', 'chr1', 1000000, 1000100):
    print(row)
```

`pileup_text` matches ordinary aligner CIGARs, including adjacent insertions
and deletions. It explicitly rejects padded (`P`) CIGARs; use `samtools
mpileup` when their text rendering is required. A deletion after a read's last
base is rendered with samtools' Q0 deletion quality rather than crashing.

## Pileup Options and Defaults

| Option | Description | Default (samtools / bcftools mpileup) | Common pitfall |
|--------|-------------|---------------------------------------|----------------|
| `-f FILE` | Reference FASTA | required by bcftools (`--no-reference` to skip) | Triggers BAQ ON by default; a missing file or contig gives reference base `N` |
| `-r REGION` | Restrict to region | | |
| `-l FILE` | BED file of regions | | |
| `-q INT` | Min mapping quality | 0 / 0 | Aligner-dependent semantics |
| `-Q INT` | Min base quality (after BAQ) | **13 / 1** | Silently drops low-quality bases, and low-quality `*` deletion slots; `-Q 0` with default overlap detection has subtle behavior |
| `--ff FLAGS` (bcftools `--ns`) | Skip reads with any of these flags | **UNMAP,SECONDARY,QCFAIL,DUP** | Pre-flagged duplicates vanish from depth; `--ff 0` keeps all mapped reads |
| `--rf FLAGS` (bcftools `--lu`) | Keep only reads with any of these flags set | none | bcftools `--nu` is not this: it skips a read missing any listed bit, i.e. keeps only reads with **all** bits |
| `-d INT` | Max depth per file | **8000 / 250** | Silently truncates; `-d 0` = no cap |
| `-B` | Disable BAQ | BAQ on with `-f` | Often correct for long reads, SV, viral, consensus; cannot combine with `-E` |
| `-A` | Count anomalous pairs (paired, not proper) | dropped | Required for amplicon (reads are by design not properly paired) |
| `-a` / `-aa` | `-a`: every position of each contig that has reads, incl. zero coverage; `-aa`: also contigs with no reads | off | Required for ARTIC SARS-CoV-2 consensus generation; `-a` = `-aa` for a single-contig reference |
| `-x` (`--disable-overlap-removal`; bcftools: `--ignore-overlaps`) | Disable mate-overlap correction | overlaps counted once | Rarely correct |
| `--max-BQ INT` (`bcftools mpileup` only) | Cap baseQ/BAQ | 60 | Not a `samtools mpileup` option; useful for ONT/HiFi (Q values inflated) |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| `[E::faidx_adjust_position] The sequence "X" was not found`, **exit status 0**, output rows with reference `N` | `-f` FASTA lacks the BAM's contig (`MN908947.3` vs `MT192765.1`, `chr1` vs `1`) | Compare `samtools view -H in.bam \| grep '^@SQ'` with `cut -f1 ref.fa.fai` before running; do not trust the exit status |
| Reference column all `N`, no error | `samtools mpileup` run without `-f` | Add `-f reference.fa` (`bcftools mpileup` refuses: "requires the --fasta-ref option") |
| `[E::mpileup] fail to parse region '22:3000-3005'` | Contig name mismatch (`chr22` vs `22`) | Use the name from the `@SQ` header |
| `[E::idx_find_and_load] Could not retrieve index file for 'in.bam'` | `-r` (or `bam.pileup()`) on an unindexed BAM | `samtools index in.bam` |
| Empty output | Region has no reads | `samtools view in.bam chr1:1000-2000 \| head`; check `ls reference.fa.fai` |
| Out of memory | High coverage region | Cap with `-d`, set well above the expected coverage (reads above the cap are dropped); or restrict with `-l targets.bed`; process contigs in parallel |
| `Failed to read from standard input: unknown file type` | Text `samtools mpileup` piped into `bcftools call` | Use `bcftools mpileup` |

Text pileup is for debugging and small regions; BCF is far cheaper to process at scale.

## Related Skills

- alignment-filtering - Filter BAM before pileup
- reference-operations - Index reference for pileup; M5 cross-check
- bam-statistics - mosdepth, depth tool selection
- variant-calling/variant-calling - Full variant calling workflows
- variant-calling/vcf-basics - VCF/BCF I/O
- variant-calling/joint-calling - Multi-sample joint calling
