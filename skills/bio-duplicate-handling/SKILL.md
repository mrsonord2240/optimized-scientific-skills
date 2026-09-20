---
name: bio-duplicate-handling
description: Mark and remove PCR/optical duplicates using samtools fixmate and markdup. Use when preparing alignments for variant calling or when duplicate reads would bias analysis.
tool_type: cli
primary_tool: samtools
license: MIT
---

## Version Compatibility

Reference examples tested with: picard 3.1+, pysam 0.22+, samtools 1.19+
Checked on: samtools 1.24, pysam 0.24.1, Picard 3.5.0, fgbio 4.1.1, umi_tools 1.1.6, samblaster 0.1.26, sambamba 1.0.1, biobambam2 2.0.185, pbmarkdup 1.2.0, mapDamage 2.2.2

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Duplicate Handling

**"Remove PCR duplicates from my BAM file"** -> Mark or remove duplicate reads using the fixmate-sort-markdup pipeline to prevent duplicate bias in variant calling.
- CLI: `samtools fixmate`, `samtools markdup` (samtools)
- Python: `pysam.fixmate()`, `pysam.markdup()` (pysam)

Mark and remove PCR/optical duplicates using samtools.

## Why Remove Duplicates?

PCR duplicates are identical copies of the same original molecule, created during library preparation. They inflate coverage, bias allele frequencies, and create false positive variant calls. Optical duplicates are flowcell-proximity artifacts: on unpatterned flowcells they arise when the imaging software splits one cluster into two adjacent calls; on patterned flowcells (NovaSeq, NovaSeq X, NextSeq 1000/2000, HiSeq X/4000) the dominant source is ExAmp (exclusion-amplification) "pad-hopping", where a library molecule re-seeds a nearby nanowell.

## When to Mark Duplicates -- and When NOT To

Standard `samtools markdup` is the right tool for some assays and actively harmful for others. The decision is assay-driven:

| Assay | Standard markdup? | Recommended approach |
|-------|------------------|----------------------|
| Germline WGS / WES (PCR or PCR-free) | YES | `samtools markdup` (PCR-free still has ~0.5% optical duplicates on patterned flowcells) |
| Somatic tumor/normal (no UMI) | YES | Same |
| Exome / target capture | YES (20-50% expected) | `samtools markdup` |
| ChIP-seq | **MARK, do not remove** | Then use peak caller's auto-dup logic (`macs3 --keep-dup auto`) |
| CUT&RUN / CUT&Tag | MARK, do not remove | Same |
| ATAC-seq | YES, BEFORE Tn5 +4/-5 shift | Then shift coords for footprinting |
| Bulk RNA-seq (no UMIs) | **NO** | Duplicates are biological at highly-expressed loci; removing them biases DE proportional to expression |
| Bulk RNA-seq (with UMIs) | NO | umi_tools dedup |
| scRNA (10x, STARsolo, drop-seq) | **NO** | umi_tools dedup with CB+UB tags, or rely on Cell Ranger UMI counts |
| ctDNA / liquid biopsy / deep panel (UMI) | NO | fgbio GroupReadsByUmi `--strategy=paired` -> CallDuplexConsensusReads (see "UMI-Aware Deduplication") |
| Twist / IDT / Roche UMI capture | NO | fgbio or Picard UmiAwareMarkDuplicatesWithMateCigar |
| Amplicon / hotspot panel (no UMI) | **NO** | Every read is a "duplicate" by coordinate; markdup erases the dataset. Use `samtools ampliconclip` instead -- see alignment-amplicon-clipping. |
| Amplicon / hotspot panel (UMI) | NO | fgbio consensus |
| Long-read native (ONT, PacBio HiFi unamplified) | NO | No PCR step; markdup is meaningless |
| PacBio HiFi amplicon | YES (rare) | `pbmarkdup` |
| Ancient DNA (aDNA) | YES + mapDamage | Run markdup, then mapDamage `--rescale` before variant calling |
| Microbiome 16S/ITS | NO | Read counts encode community structure |

If the BAM came from 10x Cell Ranger / STARsolo and `samtools markdup` produces a 50-95% duplicate rate, it is the wrong tool, not a bug.

## Tool Selection: markdup vs Picard vs UMI-aware

| Tool | Speed | Threading | Optical | UMI | Notes |
|------|-------|-----------|---------|-----|-------|
| `samtools markdup` | Fast | Yes | Yes (`-d`) | Limited (`--barcode-tag` exact-match) | Fast production choice (nf-core/sarek defaults to GATK MarkDuplicates) |
| `picard MarkDuplicates` | Slow | No | Yes | UmiAware variant (BETA, transcriptome bug) | GATK Best Practices reference |
| `biobambam2 bammarkduplicates2` | Fastest | Yes | Yes | No | Sanger / 1KGP pipelines |
| `samblaster` | Streaming, fast | No | Optional | No | Pipe directly from aligner; no name sort |
| `sambamba markdup` | Fast | Yes | Yes | No | Less actively maintained |
| `fgbio GroupReadsByUmi` + `CallMolecularConsensusReads` | Fast | Yes | n/a | **Best UMI tool** | Graph-based; supports duplex |
| `umi_tools dedup` | Slow | No | n/a | Yes (mature) | Reference for scRNA / bulk UMI |
| `pbmarkdup` | Fast | Yes | n/a | n/a | PacBio HiFi amplicons only |

Picard `UmiAwareMarkDuplicatesWithMateCigar` is BETA and has known bugs on transcriptome-aligned BAMs (silently keeps duplicates). Avoid for RNA-seq UMIs.

## Optical Distance Is Platform-Specific

`samtools markdup` default is `-d 0`, meaning **optical-duplicate detection is disabled by default**. Set explicitly per platform:

| Platform | `-d` value | Rationale |
|----------|-----------|-----------|
| HiSeq 2000/2500 (random) | 100 | Picard historic default |
| HiSeq 3000/4000/X (patterned) | 2500 | Patterned tile size larger |
| NovaSeq 6000 (patterned) | 2500 | Same as HiSeq X |
| NovaSeq X (10B) | 2500 | Patterned; same starting point as NovaSeq 6000 |
| NextSeq 1000/2000 (patterned) | 2500 | ExAmp duplicates span larger pixel distances |
| MiSeq, NextSeq 500/550 | 100 | Smaller / unpatterned |
| Element AVITI, MGI / DNBseq | Custom regex | Different read-name format -- supply via `--read-coords` |

```bash
samtools markdup -d 2500 -f stats.txt input.bam marked.bam

# Count optical (SQ) vs library/PCR (LB) duplicates. The dt:Z:SQ/LB tag is
# emitted automatically because -d is set (it is not produced by -t, which
# instead adds a 'do' tag carrying the original read's name).
samtools view -f 1024 marked.bam | grep -o 'dt:Z:[A-Z][A-Z]' | sort | uniq -c
```

Setting `-d 2500` on a HiSeq run does no harm. Forgetting `-d 2500` on NovaSeq systematically under-marks optical duplicates and overestimates library complexity.

## Multi-Library Pooled Marking

Without `--use-read-groups`, multi-library BAMs systematically over-mark: independent molecules from different libraries with the same coordinates get wrongly flagged as PCR duplicates. With `--use-read-groups`, RG tags must also match for two reads to be a duplicate (verify availability with `samtools markdup --help`):
```bash
samtools markdup --use-read-groups -d 2500 in.bam out.bam
```

samtools `--use-read-groups` keys on RG ID; Picard's library-aware behavior keys on the LB tag (allowing dedup across multiple lanes of the same library). For multi-lane single-library BAMs, Picard MarkDuplicates with `READ_NAME_REGEX` is closer to canonical.

## Duplicate Marking Workflow

**Goal:** Mark PCR/optical duplicates so they can be excluded from downstream variant calling and coverage analysis.

**Approach:** Step 0: confirm the assay against the decision table above (stop and hand off if it says NO: bulk RNA-seq, scRNA, UMI, amplicon, long-read native, 16S/ITS). Then name-sort, add mate tags with fixmate, coordinate-sort, and run markdup. The pipeline version avoids intermediate files. `examples/markdup_pipeline.sh` runs the pipeline with the assay gate, `pipefail`, a record-count check and an indexed output: `ASSAY=wgs bash examples/markdup_pipeline.sh in.bam out.bam`.

**Reference (samtools 1.19+):**
```bash
# 1. Sort by name (required for fixmate)
samtools sort -n -o namesort.bam input.bam

# 2. Add mate information with fixmate
samtools fixmate -m namesort.bam fixmate.bam

# 3. Sort by coordinate (required for markdup)
samtools sort -o coordsort.bam fixmate.bam

# 4. Mark duplicates
samtools markdup coordsort.bam marked.bam

# 5. Index result
samtools index marked.bam
```

### Pipeline Version (Optimized)
```bash
# pipefail: without it a failed first stage still exits 0 and leaves an empty marked.bam
set -euo pipefail
mkdir -p tmpdir   # collate/sort do not create it; a missing dir fails only the first stage

# collate is faster than sort -n; -u/-O between piped tools skips BGZF round-trips
samtools collate -O -u input.bam tmpdir/collate | \
    samtools fixmate -m -u - - | \
    samtools sort -u -@ 4 -T tmpdir/sort - | \
    samtools markdup -@ 4 -d 2500 --use-read-groups \
        -f markdup_stats.txt - marked.bam

samtools index marked.bam

# Sanity: markdup drops no records, and the output is not empty
test "$(samtools view -c input.bam)" -gt 0 && \
test "$(samtools view -c input.bam)" -eq "$(samtools view -c marked.bam)"
```

This is ~30% faster than `sort -n | fixmate | sort | markdup` on typical 30x WGS.

**Critical pitfall:** `samtools markdup` requires `ms` (mate score, lowercase) and `MC` (mate CIGAR) tags from `fixmate -m`. A re-sort that loses these tags (e.g. a Python round-trip) makes samtools 1.24 stop with an error (see Common Errors), not mark silently. Verify `MC:Z:` is present in the input to markdup.

**Re-marking an already-marked BAM (merged BAMs, 1000G):** add `-c` (clear previous duplicate flags and tags); without it old flags survive (111 flagged vs 101 with `-c` and Picard on a pre-marked 1000G slice).

## samtools fixmate

Adds mate information required by markdup. Must be run on name-sorted BAM.

### Basic Usage
```bash
samtools fixmate namesorted.bam fixmate.bam
```

### Add Mate Score Tag (-m)
```bash
# Required for markdup to work correctly
samtools fixmate -m namesorted.bam fixmate.bam
```

### Multi-threaded
```bash
samtools fixmate -m -@ 4 namesorted.bam fixmate.bam
```

### Remove Secondary/Unmapped
```bash
samtools fixmate -r -m namesorted.bam fixmate.bam
```

## samtools markdup

Marks or removes duplicate alignments. Requires coordinate-sorted BAM with mate tags from fixmate.

### Mark Duplicates (Keep in File)
```bash
samtools markdup input.bam marked.bam
```

### Remove Duplicates
```bash
samtools markdup -r input.bam deduped.bam
```

### Output Statistics
```bash
samtools markdup -s input.bam marked.bam 2> markdup_stats.txt
```

### Write Stats to File
```bash
samtools markdup -f stats.txt input.bam marked.bam
```

Optical distance (`-d`) is in "Optical Distance Is Platform-Specific"; threads (`-@ 4`) and `--use-read-groups` are in the pipeline. `-T PREFIX` puts temp files on a large scratch disk.

## Duplicate Statistics

### Check Duplicate Rate
```bash
samtools flagstat marked.bam
# Look for "duplicates" line
```

### Percentage Duplicates
```bash
# Primary alignments only (-F 2304 drops secondary + supplementary): the same denominator as the pysam rate below
total=$(samtools view -c -F 2304 marked.bam)
dups=$(samtools view -c -f 1024 -F 2304 marked.bam)
echo "scale=2; $dups * 100 / $total" | bc
```

A high rate means low library complexity, over-amplification or low input DNA. Over 50% on an assay that should not be marked at all means the wrong tool (see the decision table).

## pysam Python Alternative

### Full Pipeline
```python
import pysam

# Sort by name
pysam.sort('-n', '-o', 'namesort.bam', 'input.bam')

# Fixmate
pysam.fixmate('-m', 'namesort.bam', 'fixmate.bam')

# Sort by coordinate
pysam.sort('-o', 'coordsort.bam', 'fixmate.bam')

# Mark duplicates
pysam.markdup('coordsort.bam', 'marked.bam')

# Index
pysam.index('marked.bam')
```

### Check Duplicate Flag
```python
import pysam

with pysam.AlignmentFile('marked.bam', 'rb') as bam:
    total = 0
    duplicates = 0
    for read in bam:
        if read.is_secondary or read.is_supplementary:
            continue
        total += 1
        if read.is_duplicate:
            duplicates += 1

    print(f'Total: {total}')
    print(f'Duplicates: {duplicates}')
    print(f'Rate: {duplicates/total*100:.2f}%')
```

### Filter Out Duplicates
```python
import pysam

with pysam.AlignmentFile('marked.bam', 'rb') as infile:
    with pysam.AlignmentFile('nodup.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            if not read.is_duplicate:
                outfile.write(read)
```

### Production Tools, Not Hand-Rolled

For real BAMs, always use a production marker. A naive Python implementation keyed on (chrom, pos, strand) ignores 5' position correction for soft clips, ignores library/RG, treats optical = PCR, and mis-handles supplementary alignments. The result is silently wrong duplicate marks. Use `samtools markdup`, Picard, or fgbio depending on assay (see decision tables above).

## Alternative: From Aligner

Some aligners can mark duplicates directly during streaming:

### BWA-MEM2 with samblaster
```bash
# -R sets the @RG line that Picard MarkDuplicates needs (see Common Errors)
bwa-mem2 mem -R '@RG\tID:s1\tSM:s1\tLB:lib1\tPL:ILLUMINA' ref.fa R1.fq R2.fq | \
    samblaster | \
    samtools sort -o marked.bam
```

### Picard MarkDuplicates
```bash
java -jar picard.jar MarkDuplicates \
    I=input.bam \
    O=marked.bam \
    M=metrics.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
```
Picard reports READ_PAIR_DUPLICATES in pairs (samtools counts reads: 50 pairs = 100 reads). Input does not need fixmate.

### biobambam2, sambamba
```bash
# Coordinate-sorted input; no fixmate step. Metrics use Picard's column names
bammarkduplicates2 I=input.bam O=marked.bam M=metrics.txt
sambamba markdup -t 4 input.bam marked.bam
```

### Ancient DNA: mapDamage rescale after markdup
```bash
# Needs mapDamage 2.2.x (2.1.x rejects paired-end BAMs). The Bayesian rescaling step takes minutes even on small BAMs.
mapDamage -i marked.bam -r ref.fa --rescale -d mapdamage_out   # writes mapdamage_out/marked.rescaled.bam
```

### pbmarkdup (PacBio HiFi amplicons, unaligned BAM/FASTQ)
```bash
# Writes duplicate flags (0x400) into the BAM; --rmdup drops them, --dup-file keeps them in a separate file
pbmarkdup -j 4 hifi.bam marked.bam
```

## UMI-Aware Deduplication

For UMI libraries (10x scRNA, ctDNA panels, Twist/IDT/Roche UMI capture), naive markdup destroys information. Use UMI-aware tools:

### umi_tools dedup

Input must be **coordinate-sorted and indexed**.
Pass `--paired` for paired-end libraries: without it the mates are deduplicated independently and the output is silently wrong (5689 vs 2805 records on a paired-end capture BAM).

```bash
# 10x / scRNA -- group by cell barcode + UMI. Check the tags exist first: with absent CB/UB,
# --per-cell writes an EMPTY BAM and still exits 0
samtools view cellranger_possorted.bam | head -1000 | grep -c 'CB:Z:'    # must be > 0
umi_tools dedup --stdin=cellranger_possorted.bam --stdout=dedup.bam \
    --extract-umi-method=tag --umi-tag=UB --cell-tag=CB \
    --per-cell --method=directional
test "$(samtools view -c dedup.bam)" -gt 0

# Bulk UMI, paired-end (UMI in the RX tag)
samtools sort -o sorted.bam raw.bam && samtools index sorted.bam
umi_tools dedup --stdin=sorted.bam --stdout=dedup.bam --paired \
    --extract-umi-method=tag --umi-tag=RX --method=directional
```

### fgbio consensus (bulk UMI / ctDNA, best practice for low-VAF detection)

`GroupReadsByUmi` needs the mate mapping-quality (`MQ`) tag on every read (see Common Errors). `samtools fixmate -m` on name-grouped input adds it; alternatively `fgbio SetMateInformation` on queryname-sorted input. Consensus reads are written **unmapped**; re-align them before variant calling. Single-strand and duplex use different grouping strategies and are separate branches:

```bash
# If the UMI is in a separate FASTQ instead of the RX tag, annotate first and use annotated.bam below:
#   fgbio AnnotateBamWithUmis -i raw.bam -f umi.fastq -o annotated.bam
samtools sort -n -o qn.bam raw.bam
samtools fixmate -m qn.bam mated.bam        # or: fgbio SetMateInformation -i qn.bam -o mated.bam

# Single-strand molecular consensus
fgbio GroupReadsByUmi -i mated.bam -o grouped.bam --strategy=adjacency --edits=1 --raw-tag=RX
fgbio CallMolecularConsensusReads -i grouped.bam -o consensus.bam --min-reads=1

# Duplex (xGen-Prism, NEBNext duplex): needs --strategy=paired, which writes MI tags with /A /B strand
# suffixes. CallDuplexConsensusReads on adjacency-grouped reads crashes (StringIndexOutOfBoundsException).
fgbio GroupReadsByUmi -i mated.bam -o grouped_duplex.bam --strategy=paired --edits=1 --raw-tag=RX
fgbio CallDuplexConsensusReads -i grouped_duplex.bam -o duplex.bam --min-reads 1 1 0
```

### Picard UMI-aware marking
```bash
picard UmiAwareMarkDuplicatesWithMateCigar I=coordsort_fixmate.bam O=marked.bam M=metrics.txt \
    UMI_METRICS=umi_metrics.txt UMI_TAG_NAME=RX
```

`--method=directional` is the default and correct -- do not use `--method=unique`, which treats single-base UMI errors as different molecules. `samtools markdup --barcode-tag RX` (UMI/barcode handling added in samtools 1.16) does exact-match UMI grouping; adequate for IDT xGen Duplex but insufficient for single-UMI applications where 1-edit errors are common.

## Duplicate FLAG

| Flag | Value | Meaning |
|------|-------|---------|
| 0x400 | 1024 | PCR or optical duplicate |

### Filter Commands
```bash
# View only duplicates
samtools view -f 1024 marked.bam

# View non-duplicates only
samtools view -F 1024 marked.bam

# Count duplicates / non-duplicates
samtools view -c -f 1024 marked.bam
samtools view -c -F 1024 marked.bam
```

## Common Errors

Messages verbatim from samtools 1.24. Each stops the tool (exit 1); a partial output file is left behind, so delete it before re-running.

| Error | Cause | Solution |
|-------|-------|----------|
| `[bam_mating_core] ERROR: Coordinate sorted, require grouped/sorted by queryname` | fixmate input is coordinate-sorted | `samtools sort -n` (or `collate`) first |
| `samtools markdup: error, no ms score tag. Please run samtools fixmate on file first.` | fixmate was run without `-m`, or a re-sort dropped the tags | Re-run `samtools fixmate -m` |
| `samtools markdup: error, no MC tag. Please run samtools fixmate on file first.` | mate CIGAR tag lost | Re-run `samtools fixmate -m` |
| `samtools markdup: error, queryname sorted, must be sorted by coordinate.` | markdup input still name-sorted | `samtools sort` after fixmate |
| `Cannot open intermediate file "tmpdir/collate.0000.bam"` | temp dir for `collate`/`sort -T` does not exist | `mkdir -p tmpdir` |
| `Mate mapping quality (MQ) tag not present` (fgbio) | GroupReadsByUmi input has no MQ tag | `samtools fixmate -m` or `fgbio SetMateInformation` first |
| `fetch called on bamfile without index` (umi_tools) | input not sorted+indexed | `samtools sort` then `samtools index` |
| `NullPointerException ... getReadGroupId()` (Picard) | BAM has no `@RG` line | Add `@RG` (`bwa-mem2 mem -R`, `samtools addreplacerg`) |

## Lossy Operations

`samtools markdup -r` (remove duplicates) is irreversible -- the records are dropped. Default to marking, not removing; downstream tools can filter on FLAG 1024. Removing pre-emptively destroys data needed for re-running QC, library complexity estimation, or switching dedup strategies.

## Related Skills

- alignment-sorting - Sort by name/coordinate; collate vs sort -n decision
- alignment-filtering - Filter duplicates from output
- alignment-amplicon-clipping - Use ampliconclip instead of markdup for amplicon panels
- bam-statistics - Check duplicate rates with flagstat / mosdepth
- variant-calling/variant-calling - Standard variant calling expects deduped BAMs
- read-qc/quality-reports - Pre-alignment QC including UMI extraction
