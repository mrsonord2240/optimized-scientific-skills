# Selection guidance



MAPQ scales differ by aligner; the same `-q 30` filter does different things. See bio-sam-bam-basics for the full MAPQ-by-aligner table. Filtering recommendations:

| Aligner | "Drop ambiguous" | "High confidence" |
|---------|------------------|-------------------|
| BWA-MEM / BWA-MEM2 | `-q 1` | `-q 30` (or `-q 60` for unique only) |
| Bowtie2 | `-q 2` (multi-mappers get MAPQ 0 **or 1**) | `-q 23` (Bowtie2 MAPQ maxes at 42 end-to-end; 23 is a community "uniquely mapped" convention derived from analysis of Bowtie2's MAPQ scoring, not stated in the official manual) |
| **STAR** | `-q 255` | `-q 255` (STAR emits only MAPQ 0/1/3/255: 255 = unique, 3 = 2 loci, 1 = 3-4 loci, 0 = >4 loci; `-q 4` through `-q 255` are all equivalent, so `-q 60` is too) |
| HISAT2 | `-q 2` (multi-mappers get MAPQ 0 **or 1**) | `-q 60` |
| minimap2 long-read presets (`map-ont`, `map-hifi`) | `-q 1` | `-q 60` |
| minimap2 short-read (`-x sr`) | `-q 1` | `-q 30`; `-q 60` is too strict (on real Illumina reads `-q 60` kept 58.5% of primary alignments, unique MAPQ mostly 48-59) |
| pbmm2 (PacBio) | `-q 1` | `-q 60` |

For Phred-scaled aligners (BWA, minimap2), MAPQ Q maps to ~10^(-Q/10) probability of wrong mapping. For STAR, the only emitted values are 0/1/3/255 (sentinels, not probabilities).

**There is no universal "drop ambiguous" threshold.** `-q 1` removes only MAPQ 0, which is the multi-mapper value for BWA and minimap2 alone; it leaves Bowtie2 and HISAT2 multi-mappers (MAPQ 1) and STAR multi-mappers (MAPQ 1 and 3) in.
Where the aligner writes an `NH` tag (STAR, HISAT2), `-e '[NH]==1'` selects uniquely mapped alignments on any MAPQ scale. BWA, Bowtie2 and minimap2 write no `NH`, and the expression then silently returns nothing.

Checked on a 60 kb synthetic genome with planted exact and diverged repeats (BWA 0.7.19, Bowtie2 2.5.5, HISAT2 2.2.3, minimap2 2.31, STAR 2.7.11b): the "drop ambiguous" threshold above removed 400/400 exact-repeat reads for every aligner, whereas `-q 1` kept 399/400 (Bowtie2), 374/400 (HISAT2) and 80/400 (STAR). `-e '[NH]==1'` matched `-q 255` on STAR and removed 400/400 on HISAT2. On a real STAR RNA-seq BAM, MAPQ 255 coincided with `NH==1` for all 5768 records. pbmm2 26.2.99 (`--preset CCS`) on a synthetic reference with a planted 3 kb exact repeat: 5/5 reads drawn from a unique region got MAPQ 60, 4/4 reads drawn from inside the repeat got MAPQ 0.


## Assay-aware filters


**Goal:** Choose a filter that matches what the downstream caller expects. Stripping supplementary alignments breaks SV callers; requiring proper-pair drops valid spliced RNA-seq reads.

| Assay / caller | Recommended filter | Why |
|----------------|-------------------|-----|
| Germline WGS short-variant (HaplotypeCaller, DeepVariant) | `-f 2 -F 3328 -q 20` | Primary, no dup, proper pair, MAPQ>=20 |
| Somatic short-variant (Mutect2, Strelka2) | `-F 1280 -q 1` | Light pre-filter only: Mutect2 applies its own MAPQ>=20, not-secondary, not-duplicate and non-chimeric-original filters (a real Mutect2 run, GATK 4.6.2.0, removed 120 of 5642 reads at MAPQ 0 and 10 through `MappingQualityReadFilter`); supplementary reads are kept by the pre-filter because they may carry real somatic SNVs |
| Long-read short-variant (clair3, DeepVariant ONT) | `-F 3328 -q 5` | Long-read MAPQ scale is lower |
| Long-read SV (Sniffles, cuteSV) | `-F 1024` only | **Keep supplementary** -- SA tag is the SV signal |
| Short-read SV (Manta, GRIDSS, Delly, SvABA) | `-F 1024` only | Same -- supplementary required |
| ChIP-seq peak calling | `-F 1804 -q 30` | Drop dup + secondary + unmapped + mate-unmapped + QC-fail (supplementary reads are kept -- 1804 has no 2048 bit) |
| ATAC-seq | `-F 1804 -q 30 -f 2` | Same plus proper pair |
| Coverage analysis | `-F 1284 -q 1` | Mapped, no secondary, no dups; supplementary retained |
| RNA-seq quantification (STAR) | `-q 255` | Unique only (STAR sentinel) |
| RNA-seq quantification (HISAT2) | `-F 256 -q 60` | Different aligner semantics |
| RNA-seq variant (after `SplitNCigarReads`) | `-F 3328 -q 20` | Standard germline after split-N-trim |
| Panel / amplicon | After `samtools ampliconclip`; `-F 1024 -q 20` | Primer overlap makes proper-pair unreliable |
| ctDNA / cfDNA (UMI) | After fgbio consensus; do not pre-filter raw | |

The caller rationale in this table comes from caller documentation. Only the Mutect2 and HaplotypeCaller read filters were checked by running (GATK 4.6.2.0); Strelka2, DeepVariant, clair3, Sniffles, cuteSV, Manta, GRIDSS, Delly and SvABA were not run.

**Reference (samtools 1.19+):**
```bash
# Short-variant germline
samtools view -f 2 -F 3328 -q 20 -o clean.bam input.bam

# SV calling: KEEP supplementary
samtools view -F 1024 -o sv_input.bam input.bam   # NOT -F 2304, 2308, 3328 or 3332

# ChIP-seq / ATAC-seq common filter
samtools view -F 1804 -q 30 -o filtered.bam input.bam
```

**Cost of getting this wrong:** filtering out supplementary reads (`-F 2304`, `2308`, `3328` or `3332`) before SV calling produces zero SV calls -- a single-flag mistake that silently invalidates the analysis.

Flag breakdowns:
- 1280 = 256 + 1024 (secondary + duplicate)
- 1284 = 4 + 256 + 1024 (unmapped + secondary + duplicate)
- 1804 = 4 + 8 + 256 + 512 + 1024 (unmapped + mate-unmapped + secondary + QC-fail + duplicate)
- 2304 = 256 + 2048 (secondary + supplementary)
- 2308 = 4 + 256 + 2048 (unmapped + secondary + supplementary)
- 3328 = 256 + 1024 + 2048 (secondary + duplicate + supplementary)
- 3332 = 4 + 256 + 1024 + 2048 (unmapped + secondary + duplicate + supplementary)
