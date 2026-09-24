---
name: bio-alignment-amplicon-clipping
description: Trim PCR primers from aligned reads in amplicon-panel BAMs using samtools ampliconclip. Use when processing SARS-CoV-2 ARTIC, hereditary cancer panels, ctDNA hot-spot panels, or any amplicon assay where primer-derived bases would falsely confirm reference at primer footprints.
tool_type: cli
primary_tool: samtools
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: samtools 1.24, pysam 0.24.1, iVar 1.4.4 (ampliconclip needs samtools 1.11+)

Install: `conda install -c bioconda samtools pysam ivar` (pysam only for the residual-primer check, iVar only for `iVar trim`).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment Amplicon Clipping

**"Trim primer-derived bases from amplicon BAM"** -> Soft- or hard-clip the primer footprint after alignment using a primer BED, then repair fixmate/MD/NM tags.
- CLI: `samtools ampliconclip -b primers.bed input.bam -o clipped.bam` (since samtools 1.11)
- Alternative: `iVar trim` (takes a primer BED; see below)

## Why Primer Trimming After Alignment

Amplicon panels (SARS-CoV-2 ARTIC, hereditary cancer panels, ctDNA hot-spot panels, fusion panels, 16S rRNA) use designed PCR primers for enrichment. Primer-derived bases at read ends do NOT reflect biological sequence -- they reflect the primer template. Without trimming:
- False reference confirmation at primer footprint positions.
- Variant allele frequency suppressed at variants under primers (the primer sequence cannot capture the variant base).
- Strand bias artifacts (primers are typically one-strand).

Standard amplicon BAMs should NEVER be processed by `samtools markdup` -- by design every read at a primer location is a "duplicate" by coordinate. See duplicate-handling for the assay-aware decision.

## Tool Selection

| Tool | When | Notes |
|------|------|-------|
| `samtools ampliconclip` | Default for amplicon panels (since 1.11) | Soft- or hard-clip from BED; modifies CIGAR; removes MD/NM from clipped reads |
| `iVar trim` | Illumina SARS-CoV-2 / PrimalSeq route (Andersen lab) | Soft-clips by primer position + quality sliding window (`-q 20` default); expects sorted input (the index is optional in iVar 1.4.4) |
| `cutadapt` (pre-alignment) | Legacy / when alignment is downstream | Trims at FASTQ stage; less precise for amplicon |

`fgbio ClipBam` is not a primer trimmer: it clips a fixed number of bases or overlapping mate ends and takes no primer file.

## Soft-Clip vs Hard-Clip

**Goal:** Decide whether trimmed bases are kept in the BAM (reversible) or discarded (irreversible).

**Approach:** Soft-clip is the safe default; hard-clip only when archiving and disk is constrained.

| Mode | Flag | What it does | Reversible? |
|------|------|--------------|-------------|
| Soft-clip | (default) / `--soft-clip` | Bases kept in SEQ; CIGAR uses `S`; bases not aligned | Yes (CIGAR can be re-extended) |
| Hard-clip | `--hard-clip` | Bases removed from SEQ; CIGAR uses `H` | **No** (bases lost) |

Soft-clip is the recommended default. Hard-clip is irreversible -- once applied, the trimmed bases cannot be recovered for re-analysis with different primer coordinates.

## Basic ampliconclip Workflow

**Goal:** Trim primers from a coordinate-sorted, indexed amplicon BAM and produce a downstream-ready BAM.

**Approach:** Confirm the BAM is amplicon-style (read starts cluster at amplicon boundaries; `@PG` names the kit) and get the primer BED for that kit and reference build. Check contig names, run `ampliconclip` in the mode that fits the read geometry (next section), re-sort (raw output is `SO:unknown` and cannot be indexed), re-fixmate, restore MD/NM, then verify. `examples/ampliconclip_workflow.sh` runs every step and stops non-zero on any failed check; `CLIP_OPTS` sets the mode.

```bash
# 0. BAM header and BED must share a contig name (ARTIC BEDs: MN908947.3; nf-core/viralrecon
#    Illumina BAMs: MT192765.1 -- they do not pair). Empty output = mismatch, stop.
comm -12 <(samtools view -H input.bam | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2}' | sort) \
         <(awk '!/^#/{print $1}' primers.bed | sort -u)

# 1. Soft-clip primers (reversible). -f keeps the stats; "TOTAL CLIPPED: 0" means BED and BAM do not match.
samtools ampliconclip --both-ends --strand --soft-clip -f clip.stats \
    -b primers.bed input.bam -o clipped.bam
cat clip.stats

# 2. Re-sort and re-pair (CIGARs changed -- mate info needs refresh)
samtools sort -n clipped.bam | \
    samtools fixmate -m - - | \
    samtools sort -o sorted.bam -

# 3. Restore MD/NM (clipped reads lose them). Do NOT hide stderr: with a wrong reference calmd prints
#    "fail to find sequence", exits 0 and writes no MD.
samtools calmd -b sorted.bam reference.fa > clipped_final.bam
samtools index clipped_final.bam

# 4. Verify: every mapped read has MD, and no read still begins/ends inside a primer
samtools view -c -F 4 clipped_final.bam
samtools view clipped_final.bam | awk 'index($0, "\tMD:Z:"){n++} END{print n+0}'
python examples/check_primer_residual.py clipped_final.bam primers.bed --three-prime   # exit 1 if any remain
```

### Choosing the Clip Mode

Default `ampliconclip` clips only the 5' end, against any primer. Measured on the real ARTIC v5.3.2 nanopore BAM (4916 full-length amplicon reads, samtools 1.24), counting reads whose 3' end still lies inside an opposite-strand primer:

| Options | 5' end | 3' end | Reads with 3' primer left |
|---------|--------|--------|---------------------------|
| (none) | any primer | not clipped | 97.5% |
| `--strand` | primers of the read's strand | not clipped | 97.5% |
| `--both-ends` | any primer | any primer | 0% |
| `--both-ends --strand` | primers of the read's strand | primers of the opposite strand | 0% |

`--strand` matches BED column 6 to read direction and also acts together with `--both-ends`. Without it, a primer of the wrong orientation under a read clips valid sequence. Tested outcomes (synthetic single reads, primers `+` [300,325) and `-` [325,350)):

| Read | Case | none | `--strand` | `--both-ends` | `--both-ends --strand` |
|------|------|------|------------|---------------|------------------------|
| rev 251-310 | 5' end in a `+` primer | 251:50M10S (over-clipped) | untouched | 251:50M10S (over-clipped) | untouched |
| fwd 251-330 | 3' end in a `-` primer | untouched | untouched | 251:50M30S (over-clipped) | 251:75M5S |
| fwd 251-315 | 3' end in a `+` primer | untouched | untouched | 251:50M15S (over-clipped) | untouched |
| rev 331-400 | 3' end in a `-` primer | untouched | untouched | 351:20S50M (over-clipped) | untouched |

**Rule:** use `--both-ends --strand` whenever a read can reach the opposite primer: nanopore (measured) and other long reads, or any read at least as long as its amplicon. It is also safe for short reads that never reach it (synthetic paired-end panel: all 800 reads clipped to the exact planted boundary). Use `--strand` alone only when no read reaches the opposite primer; it then leaves any 3' primer in place.

Other options:
- Reads whose end matches no primer pass through unclipped and are counted as `NOT CLIPPED` in the stats (ARTIC: 92 with `--strand`, 7 with `--both-ends --strand`). `--clipped` drops them (4909 written instead of 4916), `--fail` marks them QCFAIL.
- `--tolerance N` (default 5) extends matching upstream: a read start up to N bases before the primer start still matches, while starts inside the primer already match. `--primer-counts FILE` writes reads per primer (bedgraph); `--original` adds an `OA` tag holding the pre-clip alignment; `--keep-tag` keeps the old, now wrong, NM/MD.

## Primer BED Format

```
# tab-separated, 0-based half-open like all BED; strand in column 6
chr1   100   125   primer_1_F    60   +
chr1   500   525   primer_1_R    60   -
chr1   600   625   primer_2_F    60   +
chr1   1000  1025  primer_2_R    60   -
```

Columns 1-3 give the region and column 6 the strand, which `--strand` requires: a 5-column BED fails with `Parsed 5 columns, but need at least 6`. The shipped workflow accepts tabs or spaces and ignores UCSC `track`/`browser` header lines before validating the normalized BED. ARTIC primer schemes ship pre-built BEDs (`primer.bed` from artic-network/primer-schemes) with 7 columns (chrom, start, end, name, pool, strand, sequence); ampliconclip accepts them as they are.

## SARS-CoV-2 ARTIC Comparison

| Tool | Approach | When |
|------|----------|------|
| `samtools ampliconclip` | Soft-clip from BED, post-alignment; both ends with `--both-ends --strand` | General amplicon panels; modern ARTIC workflows, nanopore included |
| `iVar trim` | Soft-clip with primer-position parsing + quality trim | nf-core/viralrecon; Illumina PrimalSeq route (Andersen lab) |

```bash
# iVar needs a coordinate-sorted BAM (the index is optional in iVar 1.4.4). -q 0 -m 1 turns the quality/length filters off so only
# primer logic acts; reads with no primer are dropped unless -e is given.
samtools index input.bam
ivar trim -i input.bam -b primers.bed -p ivar_trimmed -q 0 -m 1   # writes ivar_trimmed.bam
```

On the ARTIC v5.3.2 nanopore BAM iVar (1.4.4) wrote 4909 of 4916 reads and 4701 of them equal `ampliconclip --both-ends --strand --clipped` (name, flag, position, CIGAR), but left the 5' end inside a primer in 1.9% and the 3' end in 1.6% of reads. With its defaults (`-q 20`, `-m` half the mean read length) it wrote 0 reads on that BAM, so set `-q 0 -m 1` for nanopore data.

Note: the ARTIC network's own nanopore field-bioinformatics pipeline (`artic minion`) trims primers with its `align_trim` tool, not iVar; iVar (Grubaugh et al. 2019, Genome Biol 20:8) is the Illumina/PrimalSeq route. Modern viral consensus pipelines tend to use ampliconclip then `samtools consensus --config hiseq --ambig` (Illumina preset) for IUPAC heterozygote handling. See reference-operations for consensus generation.

## After Clipping: Required Re-Processing

Clipping changes CIGAR-derived fields and removes tags from clipped reads:

| Field | Impact | Repair |
|-------|--------|--------|
| CIGAR | New `S` or `H` operations added | Automatic from ampliconclip |
| MD:Z, NM:i | Removed from clipped reads by default (`--keep-tag` keeps the old, wrong values) | `samtools calmd -b in.bam ref.fa` recomputes both |
| Sort order | Raw output header is `SO:unknown`; `samtools index` fails with `Unsorted positions` | `samtools sort` (step 2) |
| TLEN | Template length changes when both mates clipped | `samtools fixmate -m` |
| ms, MC:Z | Mate score (lowercase per SAMtags) / mate CIGAR | `samtools fixmate -m` |

Re-calmd when a downstream tool or IGV reads MD/NM (mismatch coloring, NM-based filters). `bcftools mpileup` BAQ uses the reference, not MD: its output is identical with and without MD on the ARTIC BAM (32,150 records, checked with and without `-B`).

## Why Not Markdup

Amplicon reads at primer locations are by design coordinate-degenerate -- every read mapped to the same amplicon shares the same start/end coordinates because they all come from the same primer pair. `samtools markdup` would mark essentially every read as a duplicate and erase the dataset (synthetic paired-end panel: 710 of 800 reads). For amplicon panels:

- WITHOUT UMIs: skip dedup entirely; rely on coverage uniformity from amplicon design.
- WITH UMIs (deep panels such as Twist UMI, IDT xGen UMI, Roche AVENIO): use `fgbio GroupReadsByUmi` -> `CallMolecularConsensusReads` instead of markdup. See duplicate-handling.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Example or step 1 stats show `TOTAL CLIPPED: 0` | BAM and BED use different contig names (e.g. `MT192765.1` vs `MN908947.3`) | Use the BED for the reference the BAM was aligned to |
| `[bam_fillmd] fail to find sequence ... in the reference` and no MD tags, exit 0 | Wrong reference FASTA for calmd | Pass the FASTA the BAM was aligned to |
| `Parsed 5 columns, but need at least 6` | BED lacks the strand column | Use the 6-column format above |
| `[E::hts_idx_push] Unsorted positions` | Indexed the raw ampliconclip output | Sort first (workflow step 2) |
| High `NOT CLIPPED` share in `clip.stats` | Wrong primer-scheme version/build, or a sparse/negative-control sample | Verify scheme/build; review coverage; only raise `MAX_NOT_CLIPPED_PCT` for a documented sparse control |
| Primer bases remain at read ends (`check_primer_residual.py` exit 1, or a variant at a primer end confirms reference) | ampliconclip not run; `--strand` without `--both-ends` on reads that span the amplicon (3' primer left) | Re-run with `--both-ends --strand`; verify the clipping step ran |
| Valid sequence next to a primer is soft-clipped | `--strand` omitted, so opposite-orientation primers clip the read | Add `--strand` |
| MD/NM missing after clipping | calmd not run (ampliconclip removes them) | `samtools calmd -b clipped.bam ref.fa` |
| Markdup output shows ~100% duplicates | Amplicon BAM was processed with markdup | Restart from raw alignment; use ampliconclip; skip markdup |

## Quick Reference

| Task | Command |
|------|---------|
| Soft-clip primers (reads span amplicon, e.g. nanopore; safe default) | `samtools ampliconclip --both-ends --strand -b primers.bed in.bam -o clipped.bam` |
| Soft-clip primers (short reads that never reach the opposite primer) | `samtools ampliconclip --strand -b primers.bed in.bam -o clipped.bam` |
| Hard-clip (irreversible) | `samtools ampliconclip --both-ends --strand --hard-clip -b primers.bed in.bam -o clipped.bam` |
| Repair MD/NM after clip | `samtools calmd -b clipped.bam ref.fa > final.bam` |
| Repair mate info | `samtools sort -n - \| samtools fixmate -m - - \| samtools sort -o out.bam -` |

All ampliconclip outputs are unsorted: pipe through the repair line (or `samtools sort`) before indexing.

## Related Skills

- duplicate-handling - Why amplicon BAMs should not be markdup'd; UMI-aware alternatives
- alignment-filtering - Post-clip filtering for amplicon variant calling
- alignment-sorting - Re-sort after fixmate
- pileup-generation - depth flags for amplicon: `samtools mpileup -aa -A -d 600000 -B` (`bcftools mpileup` rejects `-aa`; use `--max-depth` and `-a FORMAT/AD,FORMAT/DP`)
- reference-operations - Consensus generation from amplicon BAMs (samtools consensus)
- read-qc/quality-reports - Pre-alignment adapter/quality trimming
