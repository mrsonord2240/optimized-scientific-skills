---
name: bio-splicing-qc
description: Assesses RNA-seq data quality for alternative splicing analysis. QC layers include experimental design audit (library prep, read length, depth, replicates), STAR cohort-style 2-pass alignment, junction saturation and plateau detection, known-vs-novel junction ratio, junction overhang and read support, splice-site strength (MaxEntScan, SpliceAI), strandedness verification, GENCODE basic vs comprehensive choice, and rRNA contamination screening. Failures in these layers silently bias PSI estimates and inflate novel-junction false positives. Use when evaluating data suitability for splicing analysis, troubleshooting low event detection, or designing sequencing experiments where AS is a primary endpoint.
tool_type: python
primary_tool: RSeQC
license: MIT
---

## Version Compatibility

Checked 2026-09-20 with: RSeQC 5.0.5, STAR 2.7.11b, samtools 1.24, pysam 0.24.1, pandas 2.3.3, maxentpy 0.0.2, Picard 3.5.0, fastq_screen 0.16.0 (minimap2 aligner), SpliceAI 1.3.1. Install (bioconda / conda-forge): `rseqc star samtools picard fastq-screen gffread maxentpy pysam pandas`; `spliceai` from PyPI. If a flag or signature differs from an older release, run `<tool> --help` or `help(module.function)` and adapt.

## Scope

Research QC of RNA-seq datasets. Splice-site scores here describe sites in a dataset; classifying a patient's variant (ACMG/ClinGen evidence, PP3/BP4) is a clinical-laboratory decision and out of scope (see `splice-variant-prediction`).

# Splicing-Specific Quality Control

Splicing analysis is more demanding than DGE on read length, depth, library prep, alignment strategy, and annotation choice. The decision sequence is: experimental design -> library prep -> alignment strategy -> annotation -> diagnostic metrics.

## Before You Run Anything

- **Gene model:** `junction_annotation.py` and `junction_saturation.py` need a **BED12** file. Convert a GTF with `gffread genes.gtf --bed | cut -f1-12 > genes.bed12` (gffread 0.12.9 appends a 13th attribute column). With a BED6 file RSeQC exits 0 and calls every junction `complete_novel`; `infer_experiment.py` accepts BED6.
- **Contig names must match the BAM** (`chr1` vs `1`). On a mismatch RSeQC exits 0 and reports zeros (`Total 0 usable reads were sampled`, flat all-zero saturation curve).
- **Rscript is needed only for RSeQC's plots.** Without it `junction_annotation.py` and `junction_saturation.py` exit 1 (`Rscript executable not found`); pass `--skip-plot` (numbers are still written to `.junction.xls` / `.junctionSaturation_plot.r`). Picard needs Java; fastq_screen needs an aligner (see rRNA section).
- **Helper script:** `examples/splicing_qc.py` wraps the RSeQC steps below with the checks above (BED12, contig overlap, no spliced reads, `--skip-plot` when Rscript is absent) and adds fragment-level junction support and MaxEntScan scoring; `python examples/test_splicing_qc.py` runs a self-contained check on a tiny planted BAM.

```bash
python examples/splicing_qc.py report sample.bam genes.bed12 sample_qc   # annotation + saturation + junction support
```

## QC Layer Taxonomy

| Layer | Tool | Section |
|-------|------|---------|
| Experimental design, library prep | Pre-sequencing review | Experimental Design Audit |
| Alignment | STAR cohort-style 2-pass | STAR 2-Pass Alignment |
| Junction discovery | RSeQC `junction_saturation`, `junction_annotation` | Junction Saturation; Novel-vs-Known |
| Junction support | pysam CIGAR parsing | Junction Read Overhang and Coverage |
| Splice site strength | MaxEntScan, SpliceAI | Splice Site Strength |
| 3' bias, mapping distribution | Picard, RSeQC `geneBody_coverage` | Picard CollectRnaSeqMetrics |
| Strand specificity | RSeQC `infer_experiment` | Strandedness Verification |
| Annotation | GENCODE basic vs comprehensive | Annotation Choice |
| Contamination | fastq_screen, samtools | rRNA Contamination Check |

All numeric pass/fail cut-offs are in **Quality Thresholds**, the design targets in **Experimental Design Audit**.

## Decision Tree by Question

| Question | Recommended QC |
|----------|-----------------|
| Will my planned RNA-seq design support AS analysis? | Pre-sequencing audit: library type, read length, depth, replicates |
| Is my data suitable for cassette exon analysis? | Junction saturation + known/novel ratio + read length |
| Why does my AS analysis call so few events? | Saturation curve, depth, library type, strandedness, 2-pass |
| Why does my AS analysis call so many novel junctions? | Annotation completeness + novel% + biology check (TDP-43, SF3B1) |
| Are my SpliceAI predictions calibrated for my tissue? | MaxEntScan + SpliceAI concordance for known sites |
| Did STAR 2-pass actually run cohort-style? | Verify the merged novel-junction file was passed to pass 2 |
| Is intron retention detectable in my data? | Library type (must be rRNA-depleted); strand-specific |

## Experimental Design Audit (Before Sequencing)

Targets are conventions, not hard limits.

| Decision | For splicing analysis | Rationale |
|----------|------------------------|-----------|
| **Library prep** | rRNA depletion (Ribo-Zero, RiboCop) | poly(A) selection loses pre-mRNA, nascent transcripts and detained introns; mandatory for IR analysis |
| **Read length** | PE 100-150 nt (PE 150 preferred) | Junction-spanning reads need >=8 nt overhang on each side |
| **Pairing** | Paired-end | Single-end loses fragment-level disambiguation of junctions |
| **Depth** | 50-100M reads/sample; >=100M for low-abundance events | 30M is the DGE-grade minimum; only a fraction of reads span junctions (24.5% on the chrX test data: 24,692 of 100,826 records) |
| **Strandedness** | Stranded library (dUTP / TruSeq stranded) | Distinguishes overlapping antisense; some tools double-count unstranded junctions |
| **Replicates** | n>=3 per condition | n=2 vs n=2 is poorly calibrated |
| **Long-intron genes (TTN, brain)** | Increase `--alignIntronMax` | Default 0 means the window-derived maximum, 2^16 x 9 = 589,824 nt (~590 kb; from `STAR --help`), which misses longer introns |

Annotation choice is in its own section below.

## STAR 2-Pass Alignment

**Goal:** Maximize novel-junction sensitivity with one junction reference shared by all samples.

**Approach:** Run STAR once per sample to discover novel junctions (pass 1), merge the novel junctions across the cohort, then re-align every sample with the merged set (pass 2). Per-sample `--twopassMode Basic` is simpler but each sample then inserts its own junctions, so junction sets differ across samples and differential calls do not replicate; use the cohort version for differential splicing. Veeneman 2016 *Bioinformatics* benchmarked per-sample two-pass (>=94% of simulated novel junctions had improved quantification); the cohort-vs-per-sample comparison is the STAR manual's multi-sample 2-pass rationale, not a number from that paper.

```bash
# Pass 1: per sample; only SJ.out.tab is needed
STAR --runMode alignReads \
    --runThreadN 8 \
    --genomeDir genome_index \
    --sjdbGTFfile gencode.v45.basic.gtf \
    --sjdbOverhang 149 \
    --readFilesIn sample_R1.fq.gz sample_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype None \
    --outFileNamePrefix pass1_${sample}_ \
    --outSJtype Standard \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 3
```

`--sjdbOverhang` = max read length - 1 (149 for 2x150). It must equal the value used at `genomeGenerate`, or STAR stops (Common Errors).

```bash
# Cohort merge. SJ.out.tab columns: 1 chr, 2 start, 3 end, 4 strand (0 undefined, 1 +, 2 -),
# 5 motif (0 non-canonical), 6 annotated (0/1), 7 unique reads, 8 multi-mapped reads, 9 max overhang.
# Keep novel (6==0), canonical-motif junctions with >=3 unique reads; STAR wants only columns 1-4.
cat pass1_*_SJ.out.tab | awk '$6 == 0 && $5 > 0 && $7 >= 3' | cut -f1-4 | sort -u > cohort_novel_SJ.tab

# Pass 2: re-align with the augmented junction set
STAR --runMode alignReads \
    --runThreadN 8 \
    --genomeDir genome_index \
    --sjdbGTFfile gencode.v45.basic.gtf \
    --sjdbFileChrStartEnd cohort_novel_SJ.tab \
    --sjdbOverhang 149 \
    --readFilesIn sample_R1.fq.gz sample_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix pass2_${sample}_ \
    --outSJtype Standard \
    --quantMode GeneCounts \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 3

samtools index pass2_${sample}_Aligned.sortedByCoord.out.bam   # pysam fetch() needs the index
```

Notes:
- `$5 > 0` alone is not a novelty filter: on the four chrX test samples the merged file had 5,717 lines but 2,184 distinct junctions, and 2,178 of those were already annotated (`sort -u` on whole lines keeps one line per differing read count). `$6 == 0` plus `cut -f1-4` leaves 6 novel junctions.
- `--outSAMstrandField intronMotif` adds the `XS` tag that leafcutter, regtools `-s XS` and Shiba need (0 of 20,000 spliced reads carried it without the flag). It also removes reads with non-canonical unannotated introns, so drop it if you need those reads.
- Pass-2 `SJ.out.tab` marks inserted novel junctions as annotated (col 6 = 1); use pass 1 to tell novel from annotated.
- Raise the unique-read threshold (or require the junction in >=N samples) if the merged file is large enough to hit the `limitSjdbInsertNsj` error.

## Junction Saturation

**Goal:** Determine whether sequencing depth is sufficient for comprehensive splicing detection.

**Approach:** RSeQC subsamples the splice events (5%, 10%, ... 100%) and counts junctions at each step. The numbers are only inside `*.junctionSaturation_plot.r` (vectors `x` percent, `y` known, `z` all, `w` novel); the helper parses them.

```bash
junction_saturation.py -i sample.bam -r genes.bed12 -o sample_junc_sat --skip-plot   # finer curve: lower -s and keep -l equal to -s (-l 2 -u 100 -s 2); -l != -s mislabels the percentages
python examples/splicing_qc.py saturation sample.bam genes.bed12 sample_junc_sat
```

**Plateau rule:** if the **known**-junction curve grows by <2% from 80% to 100% of reads, it has plateaued; still rising means more sequencing would yield more known junctions. RSeQC shuffles the events without a seed, so the curve is stochastic (three repeats on one planted BAM differed by up to 16 known junctions at the same step; growth from 80% to 100% ranged 4.0-5.3%): compare libraries by growth percentage, not by exact counts.

A curve that is flat from the first steps is a *saturated* library, not an uninformative one; a curve of all zeros means no spliced reads or a BED/BAM contig mismatch.

## Novel-vs-Known Junction Ratio

**Goal:** Detect annotation/mapping issues or biologically interesting cryptic splicing.

```bash
junction_annotation.py -i sample.bam -r genes.bed12 -o sample_junc_annot --skip-plot
```

```python
import pandas as pd

# RSeQC .junction.xls: chrom, intron_st(0-based), intron_end(1-based), read_count, annotation
junc = pd.read_csv('sample_junc_annot.junction.xls', sep='\t')
junc['annotation'] = junc['annotation'].str.strip()    # RSeQC writes ' annotated' with a leading space
by_class = junc.groupby('annotation')['read_count'].sum()
known = by_class.get('annotated', 0) / by_class.sum()
novel = (by_class.get('partial_novel', 0) + by_class.get('complete_novel', 0)) / by_class.sum()
assert abs(known + novel - 1) < 1e-9, by_class.index.tolist()
print(f'known: {known:.1%}, novel: {novel:.1%}')
```

Without the `str.strip()` the class names never match and the snippet prints `known: 0.0%, novel: 0.0%` on every library.

**Definitions.** RSeQC calls a junction `annotated` when its donor AND acceptor are each in the gene model, so an unannotated skipping junction between two annotated exons counts as known; `partial_novel` has one known end, `complete_novel` none. The snippet above is **read-weighted**, which is what the table below uses; RSeQC's printed summary is **junction-level** and reads much lower on the same BAM (planted library: 93.4% of reads known, 34.6% of junctions).

| Known fraction (reads) | Status | Interpretation |
|------------------------|--------|----------------|
| >=80% | Healthy | Comprehensive annotation, good alignment |
| 60-80% | Acceptable | Check annotation completeness or organism |
| <60% | Suspect or interesting | Mapping artifacts, contamination, OR biology |

If novel% >40%, drill down. **High novel-junction rate may be biology, not artifact:**
- **TDP-43 loss** (ALS/FTD post-mortem brain): cryptic exon de-repression in UNC13A, STMN2, ATG4B (Brown 2022 *Nature*; Klim 2019 *Nat Neurosci*)
- **SF3B1-mutant** cancer (MDS, CLL, uveal melanoma): cryptic 3'ss ~10-30nt upstream of canonical (Darman 2015 *Cell Rep*)
- **Non-model organism**: GENCODE-grade annotation unavailable; novel junctions reflect annotation gaps

## Junction Read Overhang and Coverage

**Goal:** Per-junction read support and anchor lengths, to find weakly supported junctions.

```bash
python examples/splicing_qc.py junctions sample.bam --min-overhang 8
```

`junction_stats()` in `examples/splicing_qc.py` (works on an unindexed BAM) reports, per junction `(contig, intron_start_0based, intron_end)`: `reads` (fragments passing the overhang filter), `reads_all`, and `min_overhang`. Rules it applies, each checked on planted CIGARs:
- **Overhang** is the aligned length (M/=/X) of the block immediately left and right of each N, not the total matched on that side: `30M1000N4M800N66M` has overhang 4 at both junctions (summing gives 30).
- Secondary and supplementary records are skipped; NH>1 (or MAPQ < 30 without an NH tag) is skipped; the two mates of a pair count once per junction. On real chrX pass-2 data this equals STAR's unique-read `SJ.out.tab` count on all 2,765 shared junctions with `--min-overhang 0` (371 junctions with >=10 reads, same as STAR), while per-record counting agreed on only 2,182.
- `=`/`X` CIGAR operations are handled (an `=`-blind version put a junction 50 nt off).

Junction reads with overhang <8 nt are common false positives, especially for novel sites; most callers default to >=8 nt. For very deep BAMs use STAR `SJ.out.tab` or `regtools junctions extract` instead (the helper keeps read names in memory).

## Splice Site Strength (MaxEntScan and SpliceAI)

**Goal:** Score donor and acceptor sites to flag weak / cryptic sites.

```python
import sys; sys.path.insert(0, 'examples')
from splicing_qc import score_splice_sites

# 5'ss: 9 nt (3 exon + 6 intron). 3'ss: 23 nt (20 intron + 3 exon); the intron must end in AG.
donors = ['CAGGTAAGT', 'CAGATAAGT']
acceptors = ['TTTTTTTTTTTTTTCCTTAGGAG']          # 11.58; 'T'*20 + 'CAG' has no AG and scores -7.20
s5, s3 = score_splice_sites(donors, acceptors)   # one score per input; NaN for invalid input
print(s5, s3)                                     # [10.86, 2.68] [11.58]
```

`maxentpy.maxent.score5/score3` end the process (`SystemExit: Wrong length of fa!`) on a wrong length and raise `KeyError` on N/U; lower-case is accepted. The helper validates length and A/C/G/T first and returns NaN.

| Score | Interpretation |
|-------|----------------|
| 5'ss MaxEnt > 8 | Strong donor |
| 5'ss MaxEnt 5-8 | Moderate |
| 5'ss MaxEnt < 5 | Weak / cryptic |
| 3'ss MaxEnt > 8 / < 5 | Strong / weak acceptor |

MaxEntScan (Yeo & Burge 2004 *J Comput Biol*) defines the score; the cut-offs are conventions, supported on real data: of 8,549 annotated chrX GT donors 60.7% scored >8 and 11.6% <5 (median 8.6), 95.5% of 88 decoy GT donors scored <5 (AUC 0.97 donor, 0.96 acceptor), and all 399 non-GT annotated donors scored <5.

**SpliceAI** predicts in-vivo usage from the full pre-mRNA context. Score a VCF (the assembly must match the FASTA):

```bash
spliceai -I variants.vcf -O variants.spliceai.vcf -R genome.fa -A grch38 -D 50 -M 0
# INFO: SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
```

Checked on a canonical donor G>A (PLCXD1, GRCh37): `DS_DG 0.90, DS_DL 1.00`. Reference cut-offs on the maximum delta score: >=0.2 (ClinGen SVI, Walker 2023 *Am J Hum Genet*: PP3 at supporting strength) and <=0.1 (BP4); 0.5 and 0.8 are the recommended and high-precision tiers of Jaganathan 2019 *Cell*, not ClinGen strength upgrades. Use them to prioritise sites in a dataset (Scope).

- **MaxEntScan** scores intrinsic sequence strength; **SpliceAI** scores contextual usage. High MaxEnt with low SpliceAI = intrinsically strong but contextually silenced; low MaxEnt with high SpliceAI = weak but contextually used (e.g. enhancer-driven). Report both; for variant impact see `splice-variant-prediction`.

## Picard CollectRnaSeqMetrics and Gene-Body Coverage

**Goal:** Mapping distribution (coding / UTR / intronic / intergenic / rRNA) and 5'-3' bias.

Picard needs a `refFlat` file. Build it from the BED12 (checked with Picard 3.5.0):

```bash
awk 'BEGIN{OFS="\t"} {n=$10; split($11,sz,","); split($12,st,","); s=""; e="";
     for(i=1;i<=n;i++){s=s ($2+st[i]) ","; e=e ($2+st[i]+sz[i]) ","}
     print $4,$4,$1,$6,$2,$3,$7,$8,n,s,e}' genes.bed12 > refFlat.txt

picard CollectRnaSeqMetrics \
    I=sample.bam O=sample.rna_metrics.txt REF_FLAT=refFlat.txt \
    STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
    RIBOSOMAL_INTERVALS=rRNA_intervals.interval_list      # header-bearing interval_list, not BED

geneBody_coverage.py -i sample.bam -r genes.bed12 -o sample_geneBody --skip-plot
```

Strand flags per library (foot-gun):

| Library | rMATS `--libType` | featureCounts | Picard `STRAND_SPECIFICITY` |
|---------|-------------------|---------------|------------------------------|
| Reverse-stranded (TruSeq Stranded, NEB Ultra II Directional; dUTP) | `fr-firststrand` | `-s 2` | `SECOND_READ_TRANSCRIPTION_STRAND` |
| Forward-stranded (Lexogen QuantSeq FWD, some ligation kits) | `fr-secondstrand` | `-s 1` | `FIRST_READ_TRANSCRIPTION_STRAND` |
| Unstranded | `fr-unstranded` | `-s 0` | `NONE` |

On a planted dUTP library `SECOND_READ_TRANSCRIPTION_STRAND` gave `PCT_CORRECT_STRAND_READS` 1.0 and `FIRST_READ...` 0.0. STAR has no library-strand flag; use `--outSAMstrandField intronMotif` for `XS` tags.

| Metric | Healthy (convention) | Concerning |
|--------|----------------------|------------|
| PCT_CODING_BASES | >=50% | <30% (degradation or mis-priming) |
| PCT_UTR_BASES | 20-40% | >50% (3' bias) |
| PCT_INTRONIC_BASES | <30% (poly(A)); <60% (rRNA-depleted) | >50% in poly(A): pre-mRNA contamination |
| PCT_INTERGENIC_BASES | <10% | >20%: genomic DNA contamination |
| PCT_RIBOSOMAL_BASES | see rRNA section | |
| MEDIAN_5PRIME_TO_3PRIME_BIAS | 0.7-1.3 | **<0.5 = 3' bias** (degraded RNA or poly(A) capture); **>2 = 5' bias** |
| Gene-body coverage curve | Flat | Strong 3' skew = low RIN or library mis-prep |

Checked on planted BAMs: a 3'-biased library gave `MEDIAN_5PRIME_TO_3PRIME_BIAS` 0.047 (uniform 1.009) and `geneBody_coverage` 5'/3' 0.04 vs 1.02. With 3' bias junction reads concentrate at the 3' end and miss internal junctions.

## Strandedness Verification

```bash
infer_experiment.py -i sample.bam -r genes.bed12 -s 200000     # BED6 or BED12
```

Read the two `Fraction of reads explained by ...` lines (and `failed to determine`):

| Output | Library | rMATS `--libType` |
|--------|---------|-------------------|
| ~0.5 / ~0.5 | Unstranded | `fr-unstranded` |
| >=0.9 `"1++,1--,2+-,2-+"` (PE) or `"++,--"` (SE) | Forward-stranded | `fr-secondstrand` |
| >=0.9 `"1+-,1-+,2++,2--"` (PE) or `"+-,-+"` (SE) | Reverse-stranded (dUTP / TruSeq stranded) | `fr-firststrand` |

A wrong strand setting halves usable junction reads, so verify before quantification. Leaky libraries: 20% and 40% planted leakage read 0.80 and 0.59 on the dominant string. At 0.7-0.9 report the leakage and expect that fraction of reads on the wrong strand; below 0.7 treat the library as unstranded.

## Annotation Choice

| GENCODE level | Contents | Use for |
|---------------|----------|---------|
| Basic | High-confidence canonical isoforms | Standard rMATS, leafcutter, SUPPA2; event-level AS |
| Comprehensive | All transcripts including putative/predicted | Transcript-level DTU (DRIMSeq+DEXSeq, satuRn), isoform discovery |
| RefSeq | NCBI curated | Less complete than GENCODE |
| Ensembl | Same content as GENCODE in vertebrates | Different attribute conventions |

Comprehensive captures rare isoforms but adds annotation noise and multiple-testing burden; basic can under-detect rare isoforms in DTU.

## rRNA Contamination Check

**Post-alignment (fraction of primary mapped records overlapping rRNA):**

```bash
total=$(samtools view -c -F 0x904 sample.bam)                       # -F 0x904: drop unmapped, secondary, supplementary
rrna=$(samtools view -c -F 0x904 -L rRNA_intervals.bed sample.bam)
awk -v r=$rrna -v t=$total 'BEGIN{printf "rRNA: %.1f%%\n", 100*r/t}'
```

Without `-F 0x904` the denominator includes secondary and unmapped records: a planted 22.2% library read 15.6% and passed the 20% rule. Picard `PCT_RIBOSOMAL_BASES` (0.222 on the same BAM) is the alternative.

**Pre-alignment:** `fastq_screen --aligner minimap2 --conf fastq_screen.conf --threads 8 sample_R1.fq.gz` (bowtie2, bowtie and bwa work too if installed). Setup facts, verified with fastq_screen 0.16.0:
- The default aligner is bowtie2; without it (or another aligner) the command exits 255.
- `DATABASE<TAB>rRNA<TAB>/path/to/prefix` gives the index basename. For minimap2 fastq_screen opens `prefix.fa.gz` (gzipped FASTA next to the `.mmi`); if it is missing the run exits 0 with only an `Aligner warning` in the log and reports 100% unmapped.
- Read `sample_R1_screen.txt`: the `%One_hit_one_genome` plus `%Multiple_hits_*` columns of the rRNA row are the rRNA fraction (25.00% on a planted 25% library).

| rRNA fraction | Library | Status |
|---------------|---------|--------|
| <5% | rRNA-depleted | Excellent |
| 5-20% | rRNA-depleted | Acceptable; some leakage |
| >=20% | rRNA-depleted | Failed depletion; redo |
| <5% | poly(A) | Healthy |
| >5% | poly(A) | Suggests degraded RNA or poor selection |

## Per-Tool Failure Modes

### RSeQC junction tools: silent zeros

**Trigger:** BED/BAM contig names differ, the gene model is BED6, or there are no spliced reads at MAPQ >= 30 (`-q`, default 30).

**Symptom:** exit 0 with all-zero saturation curves, `Total 0 usable reads were sampled` (`infer_experiment.py`), an empty `.junction.xls` (0 bytes, `pandas.errors.EmptyDataError`), or 0% known because every junction is `complete_novel` (BED6).

**Fix:** run the checks in **Before You Run Anything**; the helper raises a clear error for each. For `infer_experiment.py` on a BAM whose MAPQ never reaches 30 (multi-mapper-heavy aligner output), pass `-q 0`; `-q 30` is already the default.

### `infer_experiment.py`: "0 usable reads"

`Total 0 usable reads were sampled` / `Unknown data type: Mixture` means the gene model shares no contig with the BAM, the file is a GTF rather than BED, or every read fails `-q`. It is not a sample-size problem; `-s` (default 200000) only caps the sample.

### MaxEntScan: invalid sequences

See Splice Site Strength: wrong length ends the process, N/U raise `KeyError`; validate first.

## Common Errors

| Error (as printed) | Cause | Solution |
|--------------------|-------|----------|
| `EXITING because of fatal PARAMETERS error: present --sjdbOverhang=74 is not equal to the value at the genome generation step =100` | Index built with a different overhang | Use the index's value, or regenerate the index with `--sjdbOverhang` = read length - 1 |
| `Fatal LIMIT error: the number of junctions to be inserted on the fly =10552 is larger than the limitSjdbInsertNsj=...` / `SOLUTION: re-run with at least --limitSjdbInsertNsj 10552` | Annotation + merged novel junctions exceed the default 1,000,000, or a small value was set | Re-run with the value STAR prints; filter the merged file if it is unexpectedly large |
| `junction_annotation.py: error: Rscript executable not found: Rscript` (exit 1) | R not on PATH | `--skip-plot`, or install R |
| `pandas.errors.EmptyDataError: No columns to parse from file` on `.junction.xls` | BAM has no spliced reads at MAPQ >= 30 | Check the library; lower `-q` only if the aligner's MAPQ scale requires it |
| `ValueError: fetch called on bamfile without index` | pysam `fetch()` on an unindexed BAM | `samtools index`, or use `examples/splicing_qc.py` (reads with `until_eof`) |
| `SystemExit: Wrong length of fa!` / `KeyError` from `maxentpy` | 5'ss not 9 nt, 3'ss not 23 nt, or N/U present | Use `score_splice_sites` from `examples/splicing_qc.py` |

## Quality Thresholds

| Metric | Good | Acceptable | Poor | Source |
|--------|------|------------|------|--------|
| Read length (PE) | >=150 nt | 100-149 nt | <100 nt | convention |
| Sequencing depth | >=50M | 30-49M | <30M | convention |
| Junction saturation (known curve, growth 80->100%) | <2% (plateau) | Near plateau | Still rising | RSeQC convention |
| Known-junction fraction (reads) | >=80% | 60-80% | <60% (suspect or interesting) | RSeQC convention |
| Junctions with >=10 anchored reads | >=50% | 30-50% | <30% | convention; depends on depth (11% on a 50k-pair chrX sample) |
| Strandedness (dominant direction) | >=90% | 70-90% | <70% (treat as unstranded) | RSeQC convention |

Splice-site MaxEnt, rRNA and Picard cut-offs are in their own sections.

## Troubleshooting Low Event Detection

| Issue | Possible causes | Solutions |
|-------|-----------------|-----------|
| Few events called | Low depth; short reads; SE; wrong strand | Increase depth; use PE150; verify `--libType` |
| High novel junctions | Annotation gaps; mapping artifacts; biology (TDP-43, SF3B1) | Update annotation; check 2-pass; consider biology |
| Low IR detection | poly(A) library | Use rRNA depletion |
| Many weak splice sites | Cryptic splicing | Validate with MaxEnt + SpliceAI |
| PSI variance high across replicates | Library prep / RIN inconsistency | Check RIN; check 3' bias (Picard) |

## Related Skills

- splicing-quantification - PSI estimation after QC passes
- read-alignment/star-alignment - STAR 2-pass detail and parameter tuning
- read-qc/quality-reports - General sequencing QC (FastQC, MultiQC)
- read-qc/contamination-screening - rRNA / adapter / cross-species contamination
- splice-variant-prediction - SpliceAI / Pangolin for variant impact
- long-read-splicing - When short-read QC is fundamentally limiting (complex isoforms)
- differential-splicing - Downstream tool that requires QC pass

## References

- Yeo & Burge 2004 *J Comput Biol* - MaxEntScan
- Jaganathan et al 2019 *Cell* - SpliceAI
- Walker et al 2023 *Am J Hum Genet* - ClinGen SVI splicing thresholds
- Veeneman et al 2016 *Bioinformatics* - STAR 2-pass benchmark
- Brown et al 2022 *Nature* - cryptic exons in TDP-43 loss
- Klim et al 2019 *Nat Neurosci* - STMN2 cryptic splicing in ALS
- Darman et al 2015 *Cell Rep* - SF3B1 cryptic 3'ss
- Wang et al 2024 *Nat Protoc* - rMATS-turbo
- Dobin et al 2013 *Bioinformatics* - STAR aligner
