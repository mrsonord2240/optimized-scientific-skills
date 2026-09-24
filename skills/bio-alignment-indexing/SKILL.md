---
name: bio-alignment-indexing
description: Create and use BAI/CSI indices for BAM/CRAM files using samtools and pysam. Use when enabling random access to alignment files or fetching specific genomic regions.
tool_type: cli
primary_tool: samtools
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: pysam 0.22+, samtools 1.19+ (re-checked on samtools 1.24, pysam 0.24.1)

Install: `conda install -c bioconda samtools`, `pip install pysam` (pysam has no Windows wheel; use WSL/Linux).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment Indexing

Create indices for random access to alignment files using samtools and pysam.

**"Index a BAM file"** -> Create a .bai/.csi index enabling random access to genomic regions.
- CLI: `samtools index file.bam`
- Python: `pysam.index('file.bam')`

## Index Types

| Index | Extension | Max contig | Bin shift | When required |
|-------|-----------|-----------|-----------|---------------|
| BAI | `.bai` / `.bam.bai` | 2^29 bp ≈ 537 Mbp | fixed (16 kb) | Default for human, mouse, fly, fish |
| CSI | `.csi` / `.bam.csi` | 2^(min_shift + depth*3), depth sized automatically from `@SQ` lengths | `-m` sets the smallest bin (default 14) | **Required** for any contig >537 Mbp |
| CRAI | `.crai` / `.cram.crai` | chunk-based | n/a | CRAM only |
| TBI | `.tbi` | 2^29-1 | fixed | tabix VCF/BED -- same limit as BAI |

### Which Index for Which Genome

| Genome | Largest contig | Index |
|--------|---------------|-------|
| GRCh38 / GRCh37 (human) | 248 Mbp | BAI |
| GRCm39 (mouse) | 195 Mbp | BAI |
| GRCz11 (zebrafish), TAIR10 (Arabidopsis) | 78 Mbp / 30 Mbp | BAI |
| Wheat IWGSC (Triticum aestivum) | ~830 Mbp (chr3B) | **CSI** |
| Assembled pine / fir scaffolds | usually below 537 Mbp | BAI if the longest scaffold is at or below the BAI limit; check the `.fai` |
| Axolotl chromosome arms | above 537 Mbp | **CSI** (`-c`) |
| Long-read assembly with very large contigs | varies | check `cut -f2 ref.fa.fai \| sort -nr \| head -1` |

Index choice follows the **longest contig**, not total genome size. For a reference you received, check it directly:
```bash
cut -f2 ref.fa.fai | sort -nr | head -1
```
Contigs above 537 Mbp (for example wheat chromosomes or axolotl chromosome arms): use the default CSI.
```bash
samtools index -c file.bam    # min_shift 14, depth auto-sized: worked on 830 Mbp and 2.0 Gbp contigs (samtools 1.24)
```
- BAI on a contig >2^29 bp **fails loudly** and writes no `.bai`: `Region ... cannot be stored in a bai index. Try using a csi index`, exit 1.
- Do not raise `-m` for reach: depth is sized from the longest `@SQ`, and `-m` only coarsens the smallest bin (`-m 18` gave depth 4 = 2^30 on an 830 Mbp contig, not 2^33).
- BAM stores positions as int32: the header can declare a longer contig, but a read positioned beyond 2^31-1 bp (2.147 Gbp) cannot be written (`Positional data is too large for BAM format`). Split that reference before aligning; no index option fixes this.

## samtools index

### Create BAI Index
```bash
samtools index input.bam
# Creates input.bam.bai
```

### Create CSI Index
```bash
samtools index -c input.bam
# Creates input.bam.csi
```

### Specify Output Name
```bash
samtools index input.bam output.bai
```
Samtools and pysam do not auto-find a non-standard name (`Could not retrieve index file`); pass it explicitly: `samtools view -X input.bam output.bai chr1:1-1000`, or `pysam.AlignmentFile('input.bam', 'rb', index_filename='output.bai')`.

### Multi-threaded Indexing
```bash
samtools index -@ 4 input.bam
```
`-@` gave no speedup when measured (86 MB BAM: 3.1-3.2 s at `-@ 0`, 3.0 s at `-@ 4`, 3.1-3.3 s at `-@ 8`); do not expect one.

### Index CRAM
```bash
samtools index input.cram
# Creates input.cram.crai (needs no reference; reading records back does, see CRAM below)
```

## Index Requirements

Indexing requires coordinate-sorted files. The `@HD` line can be absent or wrong, so a failed index is the definitive test:
```bash
samtools view -H input.bam | grep "^@HD"    # SO:coordinate expected, but not proof
samtools index input.bam
# unsorted: [E::hts_idx_push] Unsorted positions on sequence #1: 3477 followed by 3470
#           samtools index: failed to create index (no .bai written, exit 1)

# Sort if needed, then index
samtools sort -o sorted.bam input.bam
samtools index sorted.bam
```

## Using Indices for Region Access

**Goal:** Extract reads overlapping specific genomic coordinates from an indexed BAM.

**Approach:** With the index present, `samtools view` or `pysam.fetch()` can jump directly to the relevant file offset instead of scanning the entire file.

### samtools view with Region
```bash
# Requires index file present
samtools view input.bam chr1:1000000-2000000
```

### Multiple Regions
```bash
samtools view input.bam chr1:1000-2000 chr2:3000-4000
samtools view -M input.bam chr1:1-1000000 chr1:500000-1500000   # -M merges overlapping regions
```
Without `-M`, overlapping regions return their shared reads twice (6071 vs 4555 unique in a test).

### Regions from BED File
```bash
samtools view --region-file regions.bed input.bam     # index-based (same as -M -L regions.bed)
```
Plain `samtools view -L regions.bed input.bam` is a **filter over the whole file**, not index access (0.8 s vs 0.02 s on a 600k-read BAM, same counts). It also fails on a damaged tail that region queries never read.

## CRAM

Region retrieval needs the reference; the `.crai` index alone is not enough.
```bash
samtools view -T ref.fa input.cram chr1:1000-2000        # or REF_PATH / the @SQ UR: path
```
```python
with pysam.AlignmentFile('input.cram', 'rc', reference_filename='ref.fa') as cram:
    n = cram.count('chr1', 999, 2000)
```
- Without a reachable reference `samtools view` prints no records (`Failed to populate reference`), pysam raises `OSError: truncated file`, but `samtools view -c` still counts.
- `@SQ UR:` often points at the original machine's path; pass `-T`/`reference_filename` or set `REF_PATH`.
- `REF_PATH` is an MD5 reference cache, not a directory containing `genome.fasta`: it must resolve files named from the `@SQ M5` checksum (for example, populate `REF_CACHE=cache/%2s/%2s/%s`). `-T ref.fa` is the simple route.
- pysam `get_index_statistics()` silently returns 0 for CRAM; use `samtools idxstats input.cram` or `pysam.idxstats('input.cram')`.

## pysam Python Alternative

### Create Index
```python
import pysam

pysam.index('input.bam')
# Creates input.bam.bai
```

### Create CSI Index
```python
# pysam.index passes through to samtools index; pass the -c flag for CSI.
pysam.index('-c', 'input.bam')
# Produces input.bam.csi.
```

### Fetch with Index
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    # fetch() requires index
    for read in bam.fetch('chr1', 1000000, 2000000):
        print(read.query_name)
```

### Check if Indexed (and Fresh)
```python
import os
import pysam

def index_candidates(path):
    '''Index files htslib looks for, in the order it tries them.'''
    stem = os.path.splitext(path)[0]
    if path.endswith('.cram'):
        return [path + '.crai', stem + '.crai']
    return [path + '.csi', stem + '.csi', path + '.bai', stem + '.bai']

def ensure_indexed(path):
    '''Index if missing or any index is older than the file; keep CSI if one existed.'''
    found = [i for i in index_candidates(path) if os.path.exists(i)]
    if found and all(os.path.getmtime(i) >= os.path.getmtime(path) for i in found):
        return
    csi = any(i.endswith('.csi') for i in found)
    for i in found:        # a stale .csi would still win over a fresh .bai
        os.remove(i)
    pysam.index(*(['-c'] if csi else []), path)
```
`examples/fetch_regions.py` uses this and also parses regions samtools accepts (`chr1:1,000-2,000`, `chr1:1000-`, contigs containing `:`), which `fetch(region=...)` rejects.

### Fetch Multiple Regions
```python
regions = [('chr1', 1000, 2000), ('chr1', 5000, 6000), ('chr2', 1000, 2000)]

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for chrom, start, end in regions:
        count = sum(1 for _ in bam.fetch(chrom, start, end))
        print(f'{chrom}:{start}-{end}: {count} reads')
```

### Count Reads in Region
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    count = bam.count('chr1', 1000000, 2000000)
    print(f'Reads in region: {count}')
```

### Get Reads Covering Position
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam.fetch('chr1', 1000000, 1000001):
        if read.reference_start <= 1000000 < read.reference_end:
            print(f'{read.query_name} covers position 1000000')
```

## Index File Locations

htslib (so samtools, pysam) tries, in this order, and uses the first it finds:
```
input.bam.csi   input.csi   input.bam.bai   input.bai      # BAM
input.cram.crai input.crai                                 # CRAM
```
A `.csi` therefore beats a `.bai` even when the `.bai` is newer (verified on samtools 1.24), so a stale `.csi` breaks queries after a "successful" re-index to BAI. Delete the old index first, or use `ensure_index` below. (htsjdk/GATK prefers `.bai`, the opposite order.)

## idxstats - Index Statistics

### Get Per-Chromosome Counts
```bash
samtools idxstats input.bam
```

Output format:
```
chr1    248956422    5000000    0
chr2    242193529    4500000    0
*       0            0          10000
```

Columns: reference name, length, mapped reads, unmapped reads. Counts come from the index, so a stale index reports the old totals (with an `index file is older` warning); without an index it warns `reverting to slow method` and scans the file.

### What "mapped" Actually Counts (Caveat)

The mapped column counts **every alignment record with that RNAME, including secondary AND supplementary**. For long-read minimap2 output, where a single read can produce many supplementary chimeric alignments, idxstats overcounts input reads -- typically 1.5-3x.

For unique read counts, use primary-only:
```bash
samtools view -c -F 2308 input.bam chr1   # primary mapped (2304 alone also counts placed unmapped mates)
```

Cross-check unmapped consistency (a senior sanity check):
```bash
samtools idxstats file.bam | awk '{sum+=$4} END {print sum}'   # idxstats unmapped (sum across all rows; PE orphans get a contig RNAME)
samtools view -c -f 4 -F 2304 file.bam                         # primary unmapped (should match)
```

### Sum Total Mapped Reads
```bash
samtools idxstats input.bam | awk '{sum += $3} END {print sum}'
```

### Mitochondrial Fraction
Handles `chrM` and `MT` naming and an empty BAM:
```bash
samtools idxstats input.bam | awk '$1 ~ /^(chr)?(M|MT)$/ {mt += $3} {total += $3} END {if (total) printf "MT: %.2f%%\n", mt/total*100}'
```

### pysam idxstats
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for stat in bam.get_index_statistics():
        print(f'{stat.contig}: {stat.mapped} mapped, {stat.unmapped} unmapped')
```
BAM only; for CRAM use `pysam.idxstats('input.cram')` (see CRAM above).

## FASTA Index (faidx)

Related but different - index reference FASTA for random access:

```bash
samtools faidx reference.fa
# Creates reference.fa.fai

# Fetch region from indexed FASTA
samtools faidx reference.fa chr1:1000-2000
```

### pysam FastaFile
```python
with pysam.FastaFile('reference.fa') as ref:
    seq = ref.fetch('chr1', 999, 2000)   # 0-based half-open = faidx chr1:1000-2000 (1001 bp)
    print(seq)
```

## Quick Reference

| Task | samtools | pysam |
|------|----------|-------|
| Create BAI | `samtools index file.bam` | `pysam.index('file.bam')` |
| Create CSI | `samtools index -c file.bam` | `pysam.index('-c', 'file.bam')` |
| Fetch region | `samtools view file.bam chr1:1-1000` | `bam.fetch('chr1', 0, 1000)` |
| Count in region | `samtools view -c file.bam chr1:1-1000` | `bam.count('chr1', 0, 1000)` |
| Index stats | `samtools idxstats file.bam` | `bam.get_index_statistics()` |
| Index FASTA | `samtools faidx ref.fa` | Automatic with FastaFile |

## Index Staleness

If the BAM was modified after indexing, the index points to wrong file offsets and region queries fail or return wrong reads (`The index file is older than the data file`, then `Invalid BGZF header`). Test every index that exists, not just `.bai`, and delete them all before re-indexing, because a stale `.csi` beats a fresh `.bai`. This is an mtime check only: after a restore that preserves timestamps, or if an index is truncated but newer, compare `samtools idxstats` with `samtools view -c` before trusting it.
```bash
ensure_index() {   # ensure_index file.bam|file.cram [extra samtools-index options, e.g. -@ 4]
    local f=$1; shift
    local stem=${f%.*} idx have=0 stale=0 csi=0 min_shift=14
    local -a candidates
    case $f in
        *.bam)  candidates=("$f.csi" "$stem.csi" "$f.bai" "$stem.bai") ;;
        *.cram) candidates=("$f.crai" "$stem.crai") ;;
        *) printf 'Expected a .bam or .cram file: %s\n' "$f" >&2; return 2 ;;
    esac
    for idx in "${candidates[@]}"; do
        [ -e "$idx" ] || continue
        have=1
        [ "$f" -nt "$idx" ] && stale=1
        case $idx in
            *.csi) csi=1
                    # CSI is BGZF-compressed; decompressed bytes 5-8 store min_shift.
                    min_shift=$(bgzip -cd "$idx" | od -An -j4 -N4 -tu4 | tr -d '[:space:]')
                    [ -n "$min_shift" ] || min_shift=14 ;;
        esac
    done
    if [ $have = 0 ] || [ $stale = 1 ]; then
        rm -f "${candidates[@]}"
        if [ $csi = 1 ]; then
            samtools index -c -m "$min_shift" "$@" "$f"  # preserve an existing CSI bin size
        else
            samtools index "$@" "$f"
        fi
    fi
}

shopt -s nullglob
for f in *.bam; do ensure_index "$f"; done     # an empty directory is a successful no-op
```

## Contig-Naming Sanity Check

A leading cause of "my variant calling produced empty VCFs" tickets: querying `chrM` against a BAM that uses `MT` (or `chr1` vs `1`). Always inspect contig conventions before region queries:
```bash
samtools view -H input.bam | grep '^@SQ' | head -3
# Compare with reference dict:
samtools dict ref.fa | head -3
```

UCSC convention uses `chr1`/`chrM`; Ensembl/NCBI uses `1`/`MT`. The two are not interchangeable. samtools warns and returns zero reads with exit 0; pysam raises `ValueError: invalid contig`.

## Common Errors

Messages as printed by samtools 1.24 / pysam 0.24.1.

| Error | Cause | Solution |
|-------|-------|----------|
| `Random alignment retrieval only works for indexed SAM.gz, BAM or CRAM files` + `Could not retrieve index file` (pysam: `ValueError: fetch called on bamfile without index`) | Missing index, or a non-standard index name | `samtools index file.bam`, or pass the index name (`-X` / `index_filename=`) |
| `[E::hts_idx_push] Unsorted positions on sequence #N ...` then `failed to create index` | Unsorted BAM | `samtools sort`, then index |
| `region "22:1-100" specifies an invalid region or unknown reference. Continue anyway.` (exit 0, 0 reads; pysam `invalid contig`) | Wrong contig name (`chr` vs no-`chr`, `MT` vs `chrM`) | Check `samtools view -H`; fix the name |
| `The index file is older than the data file`, then `Invalid BGZF header` / `retrieval of region failed`; zero or wrong reads on a known locus | Stale index (BAM replaced or re-sorted) | `ensure_index` above |
| `Region ... cannot be stored in a bai index. Try using a csi index`, `failed to create index` | Contig >537 Mbp with BAI | `samtools index -c file.bam` |
| `Positional data is too large for BAM format` | Contig >2^31-1 bp | Split the reference before aligning |
| `Failed to populate reference "chr1"`; pysam `OSError: truncated file` | CRAM reference unreachable | `-T ref.fa` / `reference_filename=` / `REF_PATH` |
| `fail to load index ..., reverting to slow method` (idxstats) | No index; result is still correct | Index the file |

## Related Skills

- sam-bam-basics - View and convert alignment files
- alignment-sorting - Sort BAM files (required before indexing)
- alignment-filtering - Filter by regions using index
- bam-statistics - Use idxstats for quick counts
- sequence-io/read-sequences - Index FASTA with SeqIO.index_db()
