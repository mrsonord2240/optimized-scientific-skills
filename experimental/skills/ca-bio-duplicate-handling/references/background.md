# Background

Not needed to run an analysis. This explains choices the routes and scripts already make.

## Why remove duplicates

PCR duplicates are copies of one original molecule made during library preparation. They inflate coverage, bias allele frequencies and create false-positive variant calls. Optical duplicates are flowcell-proximity artifacts: on unpatterned flowcells the imaging software splits one cluster into two adjacent calls; on patterned flowcells (NovaSeq, NovaSeq X, NextSeq 1000/2000, HiSeq X/4000) the dominant source is ExAmp "pad-hopping", where a library molecule re-seeds a nearby nanowell.

## Assay table, full

| Assay | Standard markdup? | Note |
|-------|------------------|------|
| Germline WGS / WES (PCR or PCR-free) | Yes | PCR-free still has about 0.5% optical duplicates on patterned flowcells |
| Somatic tumor/normal, no UMI | Yes | |
| Exome / target capture | Yes | 20-50% duplicates expected |
| ChIP-seq, CUT&RUN, CUT&Tag | Mark, do not remove | Peak caller handles them |
| ATAC-seq | Yes, before the Tn5 shift | |
| Ancient DNA | Yes, then mapDamage | |
| PacBio HiFi amplicon | Yes (rare) | `pbmarkdup` |
| RNA-seq, scRNA, UMI, amplicon, native long read, 16S/ITS | No | See `routes/not-markdup.md` |

## Tool comparison

| Tool | Speed | Optical | UMI | Notes |
|------|-------|---------|-----|-------|
| `samtools markdup` | Fast, threaded | `-d` | Exact-match `--barcode-tag` | nf-core/sarek defaults to GATK MarkDuplicates |
| `picard MarkDuplicates` | Slow | Yes | UmiAware variant is BETA | GATK Best Practices reference |
| `biobambam2 bammarkduplicates2` | Fast | Yes | No | Sanger / 1KGP pipelines |
| `samblaster` | Streaming | Optional | No | Pipes from the aligner; no name sort |
| `sambamba markdup` | Fast | Yes | No | Less actively maintained |
| `fgbio GroupReadsByUmi` + consensus | Fast | n/a | Best UMI tool | Graph-based; duplex supported |
| `umi_tools dedup` | Slow | n/a | Mature | Reference for scRNA and bulk UMI |

## Measured figures behind script choices

- `umi_tools dedup` without `--paired` on a paired-end capture BAM gave 5689 records against 2805 with it. Without `--random-seed=1` the count moves by about 1 between runs (5688 or 5689).
- The piped `collate | fixmate | sort | markdup` chain gives the same flagged counts as the plain `sort -n | fixmate | sort | markdup` chain and was not faster in one run (800k reads, 4 threads: 6.4 s against 4.6 s). No speed figure is claimed for 30x WGS.
- Re-marking a pre-marked 1000 Genomes slice: 111 flagged without `-c`, 101 with `-c` and with Picard.
- Setting `-d 2500` on a HiSeq run does no harm; forgetting it on NovaSeq under-marks optical duplicates and overestimates library complexity.

## Samtools versus Picard read groups

`samtools --use-read-groups` keys on the RG ID; Picard keys on the LB tag, which allows deduplication across lanes of one library. For multi-lane single-library BAMs Picard `MarkDuplicates` with `READ_NAME_REGEX` is closer to canonical.

## Open points for audit

- The source named Picard 3.1+ and samtools 1.19+ as the tested floor; the scripts were checked on samtools 1.24 and pysam 0.24.1.
- `fgbio` 4.1.1 is the version the source checked; the staged `alignment-files` environment holds a different build (see the handoff).
