---
name: bio-reference-operations
description: Generate consensus sequences and manage reference files using samtools. Use when creating consensus from alignments, indexing references, or creating sequence dictionaries.
tool_type: cli
primary_tool: samtools
license: MIT
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
- Python: iterate pileup columns and take majority base (pysam)

## samtools faidx - Index Reference FASTA

Create index for random access to reference sequences.

### Create Index
```bash
samtools faidx reference.fa
# Creates reference.fa.fai
```

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

GATK and Picard look for `<name>.dict`, where `<name>` is the FASTA file name minus `.gz` and minus its last extension: `genome.fasta` -> `genome.dict`, `ref.fa.gz` -> `ref.dict`. A `genome.fasta.dict` is ignored (GATK 4.6.2.0: `Fasta dict file .../genome.dict for reference .../genome.fasta does not exist`). `examples/prepare_reference.sh` derives the name for any of `.fa`, `.fasta`, `.fna`, `.gz`.

### With Assembly Info
```bash
samtools dict -a GRCh38 -s "Homo sapiens" reference.fa -o reference.dict
```

### Dictionary Format
```
@HD VN:1.0 SO:unsorted
@SQ SN:chr1 LN:248956422 M5:6aef897c3d6ff0c78aff06ac189178dd UR:file:reference.fa
@SQ SN:chr2 LN:242193529 M5:f98db672eb0993dcfdabafe2a882905c UR:file:reference.fa
```

The `M5:` (MD5) tag is the only definitive reference-identity check -- two references named "GRCh38" with different decoy/alt content have different M5s. CRAM enforces M5 match on read-back. See alignment-validation for BAM-vs-reference M5 cross-check.

### GRCh38 Is Not One Reference

Contig sets checked 2026-09-20 from the `.fai` / `chrom.sizes` files and NCBI's `README_analysis_sets.txt` (`ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/`).

| Reference flavor | ALT | Decoy | EBV | HLA | Use case |
|------------------|-----|-------|-----|-----|----------|
| UCSC `hg38.fa` (455 contigs) | yes (261 `_alt`) | no | no | no | UCSC browser tracks |
| NCBI `no_alt_analysis_set` | no | no | yes (`chrEBV`) | no | Aligners that are not ALT-aware |
| NCBI `no_alt_plus_hs38d1_analysis_set` | no | yes | yes | no | Same, with decoys |
| NCBI `full_analysis_set` | yes | no | yes | no | ALT-aware BWA-MEM |
| NCBI `full_plus_hs38d1_analysis_set` | yes | yes | yes | no | ALT-aware BWA-MEM with decoys |
| 1000G `GRCh38_full_analysis_set_plus_decoy_hla` = bwakit hs38DH; Broad `Homo_sapiens_assembly38.fasta` has the same 3,366 contigs | yes (261) | yes (2,385) | yes | yes (525 `HLA-*`) | GATK Best Practices, 1000G BAMs |
| T2T-CHM13 v2.0 | n/a | n/a | n/a | n/a | Distinct coordinates -- NOT interchangeable |

Mixing no-alt and ALT-aware BAMs in one cohort produces inconsistent multi-mapping behavior at HLA, KIR, and segmental-duplication regions. Standardize before joint calling.

### Contig Naming: The Silent Killer

| Convention | Source | chr1 | mitochondrion |
|-----------|--------|------|---------------|
| UCSC hg38 | UCSC Genome Browser | chr1 | chrM (16,569 bp, same sequence as Ensembl `MT`) |
| UCSC hg19 | UCSC Genome Browser | chr1 | chrM (16,571 bp, NC_001807 -- NOT the sequence of GRCh37 `MT`, 16,569 bp) |
| Ensembl (GRCh37, GRCh38) | Ensembl | 1 | MT |
| NCBI RefSeq FASTA (`GCF_*_genomic.fna`) | NCBI | NC_000001.11 (GRCh38), NC_000001.10 (GRCh37) | NC_012920.1 |
| NCBI analysis sets, 1000G GRCh38 analysis set, Broad hg38 | NCBI, 1000G, Broad | chr1 | chrM |
| 1000G phase 3 GRCh37 (`hs37d5`) | 1000G | 1 | MT (plus `GL*`, `NC_007605` EBV, `hs37d5` decoy) |

The RefSeq FASTA has no `chr1`; the `_assembly_report.txt` beside it maps every name (tab-separated, CRLF line ends): column 1 Sequence-Name (`1`, `X`, `MT` = Ensembl names for the chromosomes only; scaffolds get GRC names there, Ensembl uses their GenBank accession, e.g. `KI270706.1`), 5 GenBank (`CM000663.2`), 7 RefSeq (`NC_000001.11`), 10 UCSC-style-name (`chr1`). Use it as the rename map for `_alt`, `_random` and `chrUn` contigs; stripping `chr` only covers chr1-22, X, Y, M.

A BAM with `@SQ SN:chr1` cannot be analyzed against a `1`-named reference (and vice versa). Detect:
```bash
samtools view -H sample.bam | grep '^@SQ' | head -3
samtools dict ref.fa | head -3
```

### Rename Contigs Without Re-aligning

Renaming is a header-only change: BAM records store a contig index, so `samtools reheader` gives records identical to the input except the name. Only do it when the sequences are the same (compare `LN`, and `M5` where the BAM header has it -- a UCSC hg19 `chrM` is not GRCh37 `MT`).
```bash
# map.tsv: old<TAB>new, one contig per line (e.g. chr22<TAB>22, chrM<TAB>MT), or UCSC -> RefSeq from the assembly report:
grep -v '^#' GCF_000001405.40_GRCh38.p14_assembly_report.txt | tr -d '\r' | awk -F'\t' '$10!="na"{print $10 "\t" $7}' > map.tsv
samtools view -H sample.bam | awk -F'\t' -v OFS='\t' 'NR==FNR{m[$1]=$2; next}
    /^@SQ/{for(i=2;i<=NF;i++) if($i~/^SN:/){n=substr($i,4); if(n in m) $i="SN:" m[n]}} {print}' map.tsv - > renamed.hdr
samtools reheader renamed.hdr sample.bam > renamed.bam && samtools index renamed.bam

# UCSC -> Ensembl for the primary chromosomes only (hg38/GRCh38, not hg19):
samtools view -H sample.bam | sed -e 's/^\(@SQ\tSN:\)chrM/\1MT/' -e 's/^\(@SQ\tSN:\)chr/\1/' > renamed.hdr

# Check: contig names and lengths now equal the reference's
diff <(samtools view -H renamed.bam | awk '/^@SQ/{print $2, $3}') <(awk '{print "SN:"$1, "LN:"$2}' ref.fa.fai) && echo OK
```
`SA:Z:`, `XA:Z:` and `OA:Z:` tags hold contig names as text and keep the old names; drop them with `samtools view -b -x SA -x XA renamed.bam` if a downstream tool reads them. `samtools reheader` also accepts a CRAM (checked on 1.24; it warns `Failed to populate reference` when no `REF_PATH`/cache holds the renamed reference, and the decode with `-T` is identical). For VCF use `bcftools annotate --rename-chrs map.tsv`.

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

### IUPAC Ambiguity for Heterozygotes
```bash
# Default Bayesian mode. --ambig is REQUIRED for IUPAC codes (R, Y, S, W, K, M, B, D, H, V);
# without it an ambiguous column is N, even when one base is 80% of the reads
samtools consensus --ambig input.bam -o consensus.fa

# Tune Bayesian het calling with --het-scale (< 1 fewer IUPAC calls, > 1 more)
samtools consensus --ambig --het-scale 0.1 input.bam -o consensus.fa

# Fixed fractional thresholds exist only in -m simple
samtools consensus -m simple --ambig --het-fract 0.2 --call-fract 0.5 input.bam -o consensus.fa
```

`--het-fract` and `--call-fract` are ignored in the default Bayesian mode (byte-identical output for 0.05, 0.9 and 0.2/0.5) and act only with `-m simple`:
- `--het-fract F`: minimum ratio of the second-most to the most common base for an IUPAC call (needs `--ambig`). Always pass it explicitly: with the flag omitted, a column with 15% minor allele stays a plain base although the help prints a default of 0.15.
- `--call-fract F`: fraction of reads that must agree on the top base, otherwise `N` (default 0.75).

`--show-ins` / `--show-del` control insertion / deletion display, not ambiguity.

### Platform-Aware Consensus
```bash
# The default Bayesian algorithm needs no --config; platform-specific profiles (samtools 1.17+; list them with `samtools help consensus`)
samtools consensus --config hifi       input.bam -o consensus.fa   # PacBio HiFi
samtools consensus --config r10.4_sup  input.bam -o consensus.fa   # ONT R10.4+ (r10.4_dup for duplex)
samtools consensus --config ultima     input.bam -o consensus.fa   # Ultima Genomics
samtools consensus --config hiseq      input.bam -o consensus.fa   # Illumina

# Report ref base where consensus unavailable (low coverage; -T added in samtools 1.22; bases keep the FASTA's case)
samtools consensus -T ref.fa input.bam -o consensus.fa
```

### samtools consensus vs bcftools consensus

Different operations -- conflating them produces nonsense:

| Tool | Input | Output | Use case |
|------|-------|--------|----------|
| `samtools consensus` | BAM | Consensus FASTA derived from reads (Bayesian) | Viral, de novo / amplicon, low-coverage species |
| `bcftools consensus` | reference + VCF | Reference with VCF variants applied | Apply called variants (haplotype reconstruction, custom ref for re-mapping) |

For viral consensus from BAM:
```bash
# Modern: samtools consensus
# (--show-del yes would write '*' into the FASTA, so the default no is kept)
samtools consensus --config hiseq -d 10 --ambig -a input.bam -o consensus.fa

# Apply called variants to reference (different question)
bcftools consensus -f reference.fa variants.vcf.gz -o sample_consensus.fa
bcftools consensus -f reference.fa -H 1 phased.vcf.gz -o haplotype1.fa   # phased haplotype 1
```

With a genotyped (FORMAT/GT) VCF and no `-H`, bcftools 1.24 writes heterozygous SNPs as IUPAC codes; use `-H 1` / `-H 2` for one haplotype, or `-H A` (or `-s -`) to apply every ALT allele.

`samtools consensus` is not iterative and is not an assembly-polishing tool.

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
        end = min(end, ref.get_reference_length(chrom))
        if start >= end:
            print(f'skip {chrom}:{start + 1}-{end}: outside the contig')
            continue
        seq = ref.fetch(chrom, start, end)
        print(f'>{chrom}:{start + 1}-{end}')
        for i in range(0, len(seq), 60):
            print(seq[i:i + 60])
```

### Generate Simple Consensus
```python
import pysam
from collections import Counter

def consensus_at_position(bam, chrom, pos):
    bases = Counter()
    for pileup in bam.pileup(chrom, pos, pos + 1, truncate=True):
        if pileup.pos == pos:
            for read in pileup.pileups:
                if not read.is_del and not read.is_refskip:
                    bases[read.alignment.query_sequence[read.query_position]] += 1
    if bases:
        return bases.most_common(1)[0][0]
    return 'N'

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    consensus = consensus_at_position(bam, 'chr1', 1000000)   # 0-based position
    print(f'Consensus at chr1:{1000000 + 1} = {consensus}')
```

### Build Consensus Sequence (Pedagogical Only)

The Python majority-vote consensus below is illustrative, NOT production. `samtools consensus` is Bayesian, quality-aware, and platform-aware; majority vote weights every base equally and produces wrong calls on low-coverage / low-quality regions. Use for teaching pileup iteration mechanics; use `samtools consensus` for any real consensus.

`pileup()` filters before the vote (pysam 0.24.1 defaults): bases with quality < 13, unmapped / secondary / QC-fail / duplicate reads, and orphan reads are dropped, and overlapping mates are counted once. `max_depth` defaults to 8000 (deeper columns are subsampled), so raise it.

`build_consensus` returns exactly `end - start` characters, so index `i` is reference position `start + i`. `pileup()` skips uncovered columns; the function starts from all-`N` and fills only the columns it sees. Building the string by appending per pileup column shifts everything after the first coverage gap (581 false differences vs 1 true on the real chr22 slice).

```python
import pysam
from collections import Counter

def build_consensus(bam_path, chrom, start, end, min_depth=3):
    """Majority vote over [start, end), 0-based half-open; N where depth < min_depth."""
    consensus = ['N'] * (end - start)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup in bam.pileup(chrom, start, end, truncate=True, max_depth=1_000_000):
            bases = Counter()
            for read in pileup.pileups:
                if not read.is_del and not read.is_refskip:
                    base = read.alignment.query_sequence[read.query_position]
                    bases[base.upper()] += 1

            if sum(bases.values()) >= min_depth:
                consensus[pileup.reference_pos - start] = bases.most_common(1)[0][0]

    return ''.join(consensus)
```

### Compare Consensus to Reference (Python)
```python
def compare_to_ref(bam_path, ref_path, chrom, start, end, min_depth=3):
    """[(1-based position, ref base, consensus base)] for called bases that differ from the reference."""
    consensus = build_consensus(bam_path, chrom, start, end, min_depth)
    with pysam.FastaFile(ref_path) as ref:
        reference = ref.fetch(chrom, start, end).upper()   # soft-masked FASTA is lowercase
    return [(start + i + 1, r, c)
            for i, (c, r) in enumerate(zip(consensus, reference))
            if c != 'N' and c != r]
```
Ties (50/50 columns) go to the first base counted; `samtools consensus` calls them `N` (or an IUPAC code with `--ambig`) and weights bases by quality, so expect a few different calls at het columns and at shallow, low-quality columns (chr22 slice: 1 difference by majority vote and by `-m simple --call-fract 0.5 --min-BQ 13`, 2 by the default Bayesian mode).

### Header Dict for Writing a BAM (not a .dict file)
`pysam.AlignmentFile(..., 'wb', header=header)` takes this dict. It has no `M5`, so it is not a sequence dictionary: use `samtools dict` for that.
```python
import pysam

def create_dict_header(fasta_path):
    header = {'HD': {'VN': '1.6', 'SO': 'unsorted'}, 'SQ': []}

    with pysam.FastaFile(fasta_path) as ref:
        for name in ref.references:
            length = ref.get_reference_length(name)
            header['SQ'].append({'SN': name, 'LN': length})

    return header

header = create_dict_header('reference.fa')
for sq in header['SQ'][:5]:
    print(f'{sq["SN"]}: {sq["LN"]:,} bp')
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
For a per-position list of differences, use `compare_to_ref` above.

## Related Skills

- sam-bam-basics - CRAM reference resolution (REF_PATH, REF_CACHE)
- alignment-indexing - faidx for reference access
- alignment-validation - BAM-vs-reference M5 cross-validation
- pileup-generation - Pileup for consensus building
- variant-calling/vcf-basics - VCF I/O for `bcftools consensus`
- variant-calling/consensus-sequences - Consensus from VCF (different operation)
- read-alignment/bwa-alignment - BWA index preparation
- sequence-io/read-sequences - Parse FASTA with Biopython
