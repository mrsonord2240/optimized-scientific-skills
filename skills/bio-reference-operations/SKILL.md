---
name: bio-reference-operations
description: Generate consensus sequences and manage reference files using samtools. Use when creating consensus from alignments, indexing references, creating sequence dictionaries, extracting regions, renaming or matching contig names (chr22 vs 22, GRCh38 flavours), or resolving CRAM references.
tool_type: cli
primary_tool: samtools
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: GATK 4.5+, bcftools 1.19+, pysam 0.22+, samtools 1.19+
Checked 2026-09-20 on samtools 1.24, bcftools 1.24, pysam 0.24.1, GATK 4.6.2.0, Picard 3.5.0.
Install: `conda install -c bioconda samtools bcftools` and `pip install pysam`.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Reference Operations

Generate consensus sequences and manage reference files using samtools.

**"Prepare a reference genome"** -> Index the FASTA and create a sequence dictionary for downstream tools.
- CLI: `samtools faidx ref.fa` + `samtools dict ref.fa -o ref.dict`
- Python: `pysam.FastaFile('ref.fa')` (auto-uses .fai index)

**"Build a consensus from BAM"** -> Derive the most-supported base at each position from aligned reads.
- CLI: `samtools consensus input.bam -o consensus.fa`
- Python: iterate pileup columns and take majority base (pysam), `references/python-consensus.md`
- IUPAC codes, `--het-fract`, platform `--config`, `-T`, viral consensus, `bcftools consensus`: `references/consensus-modes.md`

**"BAM says `chr22`, reference says `22`" / "which GRCh38 is this?"** -> rename contigs with `samtools reheader` (no re-alignment), or identify the reference flavour: `references/contig-naming.md`.

## samtools faidx - Index Reference FASTA

Create index for random access to reference sequences.

### Create Index
```bash
samtools faidx reference.fa
# Creates reference.fa.fai
```
A `.gz` FASTA must be bgzip-compressed; plain gzip (Ensembl/NCBI downloads) fails with `Cannot index files compressed with gzip, please use bgzip`. Recompress with `gunzip -c ref.fa.gz | bgzip > ref.bgz.fa.gz`.

### Fetch Region from Reference
```bash
samtools faidx reference.fa chr1:1000-2000
```

### Fetch Multiple Regions
```bash
samtools faidx reference.fa chr1:1000-2000 chr2:3000-4000
```

### Fetch Entire Chromosome
```bash
samtools faidx reference.fa chr1
```

### Output to File
```bash
samtools faidx reference.fa chr1:1000-2000 > region.fa
```

### Reverse Complement
```bash
samtools faidx -i reference.fa chr1:1000-2000                    # header becomes ">chr1:1000-2000/rc"
samtools faidx -i --mark-strand no reference.fa chr1:1000-2000   # keep the plain header
```

### faidx Errors (samtools 1.24)
`faidx` builds the `.fai` itself the first time a region is requested; there is no need to index first.

| Situation | Message | Exit |
|-----------|---------|------|
| Contig name not in the FASTA (`22` vs `chr22`) | `[faidx] Failed to fetch sequence in 22:1-100` | 1 |
| FASTA file missing | `Failed to open the file ... Could not load fai index` | 1 |
| Region entirely beyond the contig end | `[faidx] Zero length sequence` and an empty record | **0** |

Check names with `cut -f1 reference.fa.fai` (names are the first word of each `>` line), and check coordinates against the `.fai` lengths before trusting a record that has a header but no bases.

### FAI File Format
```
chr1    248956422    6    60    61
chr2    242193529    253105708    60    61
```
Columns: name, length, offset, line bases, line width

## samtools dict - Create Sequence Dictionary

Create SAM header dictionary for reference (used by GATK, Picard).

### Create Dictionary
```bash
samtools dict reference.fa -o reference.dict
```

GATK and Picard look for `<name>.dict`, where `<name>` is the FASTA file name minus `.gz` and minus its last extension: `genome.fasta` -> `genome.dict`, `ref.fa.gz` -> `ref.dict`. GATK ignores `genome.fasta.dict` (GATK 4.6.2.0 HaplotypeCaller: `Fasta dict file .../genome.dict for reference .../genome.fasta does not exist`); Picard 3.5.0 accepts both names, so `<name>.dict` is the safe one for both. `examples/prepare_reference.sh` derives the name for any of `.fa`, `.fasta`, `.fna`, `.gz`.

### With Assembly Info
```bash
samtools dict -a GRCh38 -s "Homo sapiens" reference.fa -o reference.dict
```

### Dictionary Format
```
@HD VN:1.0 SO:unsorted
@SQ SN:chr1 LN:248956422 M5:6aef897c3d6ff0c78aff06ac189178dd UR:file:///abs/path/to/reference.fa
@SQ SN:chr2 LN:242193529 M5:f98db672eb0993dcfdabafe2a882905c UR:file:///abs/path/to/reference.fa
```
`UR` is the absolute `file:///` path of the FASTA (`-u` overrides it). `examples/toy.fa` and `examples/toy.expected.dict` are a two-contig check: `samtools dict -u toy.fa toy.fa | diff - toy.expected.dict` prints nothing.

The `M5:` (MD5) tag is the only definitive reference-identity check -- two references named "GRCh38" with different decoy/alt content have different M5s. CRAM enforces M5 match on read-back. See alignment-validation for BAM-vs-reference M5 cross-check.

## samtools consensus - Generate Consensus

Create consensus sequence from alignments.

### Basic Consensus
```bash
samtools consensus input.bam -o consensus.fa
```

### From Specific Region
```bash
samtools consensus -r chr1:1000-2000 input.bam -o region_consensus.fa
```

### Output Formats
```bash
# FASTA (default)
samtools consensus -f fasta input.bam -o consensus.fa

# FASTQ (includes quality); wrapped at 70 columns unless -l 0 (then 4 lines per record)
samtools consensus -f fastq -l 0 input.bam -o consensus.fq
```

### Quality Options
```bash
# Minimum depth to call base (positions below it become N)
samtools consensus -d 5 input.bam -o consensus.fa

# Pad the start/end of each contig with N up to the reference ends (header LN)
samtools consensus -a input.bam -o consensus.fa
```

Columns inside the covered span are always emitted (N where nothing is called), so `-a` changes only the ends: `-a` pads contigs that have reads, `-aa` also emits contigs with no reads (all N). Calls inside the span are identical with or without it. For exactly one character per reference position (coordinates line up with the reference) use `-a --show-del yes --show-ins no`; deleted columns become `*`.

Option list: run `samtools help consensus` (exit 0; `samtools consensus --help` prints the usage after an "unrecognized option" error and exits 1) or `man samtools-consensus`.

## pysam Python Alternative

### Fetch from Indexed FASTA
```python
import pysam

with pysam.FastaFile('reference.fa') as ref:
    seq = ref.fetch('chr1', 999, 2000)  # 0-based
    print(seq)
```

### Get Reference Lengths
```python
with pysam.FastaFile('reference.fa') as ref:
    for name in ref.references:
        length = ref.get_reference_length(name)
        print(f'{name}: {length:,} bp')
```

### Fetch All Chromosomes
```python
with pysam.FastaFile('reference.fa') as ref:
    for chrom in ref.references:
        seq = ref.fetch(chrom)
        print(f'>{chrom}')
        print(seq[:100] + '...')
```

### Fetch Several Regions (0-based)
pysam coordinates are 0-based half-open; samtools regions are 1-based inclusive (`fetch('chr1', 999, 2000)` == `faidx chr1:1000-2000`). `fetch` silently truncates at the contig end and returns `''` past it, so bound the request:
```python
regions = [('chr1', 0, 10000), ('chr2', 5000, 15000)]
with pysam.FastaFile('reference.fa') as ref:
    for chrom, start, end in regions:
        length = ref.get_reference_length(chrom)
        if start >= length:
            print(f'skip {chrom}:{start + 1}-{end}: outside the contig ({length} bp)')
            continue
        end = min(end, length)
        seq = ref.fetch(chrom, start, end)
        print(f'>{chrom}:{start + 1}-{end}')
        for i in range(0, len(seq), 60):
            print(seq[i:i + 60])
```

## Reference Preparation Workflow

**Goal:** Set up a reference genome with all indices needed by common analysis tools.

**Approach:** Create FASTA index (.fai), sequence dictionary (.dict), and aligner-specific indices in sequence.

### Prepare Reference for Analysis
Index (`samtools faidx`) and dictionary (`samtools dict`, named `<name>.dict`) as above; for any of `.fa`, `.fasta`, `.fna`, `.gz`, plus chromosome sizes: `bash examples/prepare_reference.sh reference.fa`. Then, for CRAM:
```bash
# Pre-populate CRAM REF_CACHE (for offline HPC nodes); the cache dir must exist and be yours
REF_CACHE_DIR=$HOME/ref_cache
mkdir -p "$REF_CACHE_DIR"
seq_cache_populate.pl -root "$REF_CACHE_DIR" reference.fa

# Decode a CRAM with no -T and no network
REF_PATH="$REF_CACHE_DIR/%2s/%2s/%s" samtools view input.cram | head
```

Without `REF_PATH`/`REF_CACHE` a CRAM decode needs the reference: `samtools view -T reference.fa input.cram`.

For aligner-specific indices (BWA, Bowtie2, STAR, minimap2, Salmon), see read-alignment.

### Check Reference Setup
```bash
REF=reference.fa
NAME=$(basename "${REF%.gz}"); DICT="$(dirname "$REF")/${NAME%.*}.dict"   # genome.fasta -> genome.dict
[ -s "${REF}.fai" ] && echo "FAI: OK" || echo "FAI: MISSING"
[ -s "$DICT" ] && echo "DICT: OK ($DICT)" || echo "DICT: MISSING ($DICT)"
samtools faidx "$REF" "$(head -1 "${REF}.fai" | cut -f1):1-100" > /dev/null && echo "Fetch: OK" || echo "Fetch: FAILED"
```

## Common Operations

### Extract Chromosome
```bash
samtools faidx reference.fa chr1 > chr1.fa
samtools faidx chr1.fa  # Index the subset
```

### Get Chromosome Sizes
```bash
cut -f1,2 reference.fa.fai > chrom.sizes
```

### Subset Reference
```bash
# faidx exits 1 on a missing contig (Ensembl names `1` vs `chr1`), so && stops the chain
samtools faidx reference.fa chr1 chr2 chr3 > subset.fa &&
    samtools faidx subset.fa &&
    samtools dict subset.fa -o subset.dict
```

### Compare Consensus to Reference
```bash
# Generate consensus
samtools consensus input.bam -o consensus.fa

# Align consensus back to reference
minimap2 -a reference.fa consensus.fa > comparison.sam
```
For a per-position list of differences, use `compare_to_ref` in `references/python-consensus.md`.

## Reference Files

| File | Read when |
|------|-----------|
| `references/contig-naming.md` | BAM and reference disagree on contig names, "which GRCh38 is this", renaming contigs with `samtools reheader`, UCSC / Ensembl / RefSeq name maps |
| `references/consensus-modes.md` | IUPAC codes (`--ambig`, `--het-fract`, `--call-fract`, `--het-scale`), platform `--config` profiles, `-T`, viral consensus, `samtools consensus` vs `bcftools consensus` |
| `references/python-consensus.md` | pysam majority-vote consensus, `compare_to_ref`, the header dict for writing a BAM |
| `scripts/rename_contigs.sh` | Rename BAM contigs from a `map.tsv` (`in.bam map.tsv out.bam [ref.fa]`), see `references/contig-naming.md` |
| `scripts/pysam_consensus.py` | Majority-vote consensus / differences from the reference for a BAM window (`consensus` and `compare` subcommands) |

## Related Skills

- sam-bam-basics - CRAM reference resolution (REF_PATH, REF_CACHE)
- alignment-indexing - faidx for reference access
- alignment-validation - BAM-vs-reference M5 cross-validation
- pileup-generation - Pileup for consensus building
- variant-calling/vcf-basics - VCF I/O for `bcftools consensus`
- variant-calling/consensus-sequences - Consensus from VCF (different operation)
- read-alignment/bwa-alignment - BWA index preparation
- sequence-io/read-sequences - Parse FASTA with Biopython
