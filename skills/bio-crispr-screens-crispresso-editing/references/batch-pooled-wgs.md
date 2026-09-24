# Batch, pooled-amplicon and WGS modes (CRISPResso2)

Moved from SKILL.md. Read when the design is many samples (CRISPRessoBatch), many amplicons in one library (CRISPRessoPooled) or an off-target survey from a BAM (CRISPRessoWGS).

## Batch Mode (Multi-Sample, Same Amplicon)

**Goal:** Process tens to hundreds of samples with same amplicon design (e.g., a timecourse, dose response, or replicate panel).

**Approach:** Provide a tab-separated batch settings file with per-sample parameters; CRISPRessoBatch runs all in parallel.

```bash
# batch_settings.txt (tab-separated, headers required)
# name    fastq_r1                fastq_r2                amplicon_seq    guide_seq
# t0      t0_R1.fq.gz             t0_R2.fq.gz             ACGT...         GUIDE
# t6      t6_R1.fq.gz             t6_R2.fq.gz             ACGT...         GUIDE
# t12     t12_R1.fq.gz            t12_R2.fq.gz            ACGT...         GUIDE
# t24     t24_R1.fq.gz            t24_R2.fq.gz            ACGT...         GUIDE

CRISPRessoBatch \
    --batch_settings batch_settings.txt \
    --batch_output_folder batch_run \
    --skip_failed \
    --n_processes 8

# Outputs:
#   batch_run/CRISPRessoBatch_on_<batch file name>/CRISPRessoBatch_RUNNING_LOG.txt
#   batch_run/CRISPRessoBatch_on_<batch file name>/CRISPRessoBatch_quantification_of_editing_frequency.txt  (aggregated)
#   batch_run/CRISPRessoBatch_on_<batch file name>/CRISPResso_on_<name>/ for each sample
```

## Pooled-Amplicon Mode

**Goal:** Process multi-amplicon sequencing libraries (e.g., arrayed validation pools).

**Approach:** Provide an amplicon table with one row per target; CRISPRessoPooled de-multiplexes reads to the correct amplicon.

```bash
# amplicons.txt (tab-separated; header may vary by CRISPResso2 version)
# amplicon_name  amplicon_seq    guide_seq
# BRCA1_exon3    ACGT...         GUIDE1
# TP53_exon7     ACGT...         GUIDE2
# KRAS_codon12   ACGT...         GUIDE3

# --min_reads_to_use_region 100: see the note below (the default of 1000 skips small amplicons)
CRISPRessoPooled \
    --fastq_r1 pooled_R1.fastq.gz \
    --fastq_r2 pooled_R2.fastq.gz \
    --amplicons_file amplicons.txt \
    --output_folder pooled_run \
    --min_reads_to_use_region 100 \
    --n_processes 8

# Outputs:
#   pooled_run/CRISPRessoPooled_on_<fastq name>/SAMPLES_QUANTIFICATION_SUMMARY.txt
#   pooled_run/CRISPRessoPooled_on_<fastq name>/CRISPResso_on_<amplicon>/ for each amplicon
```

**`--min_reads_to_use_region` defaults to 1000.** Any amplicon with fewer aligned reads than this is silently
skipped -- CRISPRessoPooled still exits 0 and writes a complete-looking `SAMPLES_QUANTIFICATION_SUMMARY.txt`,
but every field for that amplicon is `NA` (confirmed: a real 2-amplicon, ~250-reads/amplicon pilot pool --
exactly the "arrayed validation pool" use case this mode is for -- returns all-NA at the default). Set
`--min_reads_to_use_region` below the expected per-amplicon read depth for pilot/validation-scale pools.
Always check `SAMPLES_QUANTIFICATION_SUMMARY.txt` for `NA` rows before trusting the output.

## WGS Off-Target Mode

**Goal:** Quantify off-target editing from whole-genome sequencing.

**Approach:** Provide BAM file + reference + BED file of suspected off-target sites; CRISPResso extracts reads from each region and quantifies edits.

```bash
CRISPRessoWGS \
    --bam_file aligned.bam \
    --reference_file genome.fa \
    --region_file off_targets.txt \
    --output_folder wgs_run \
    --n_processes 8

# Outputs:
#   wgs_run/CRISPRessoWGS_on_<bam name>/SAMPLES_QUANTIFICATION_SUMMARY.txt
#   wgs_run/CRISPRessoWGS_on_<bam name>/CRISPResso_on_<region>/ for each analysed region
```

The flag names are exact: `--bam` alone is rejected (`ambiguous option: --bam could match --bam_output, --bam_file`).
`--region_file` is a tab-separated `chr  start  end  name` table with no header (e.g. `chr9  962  1198  HEK3`).
**`--min_reads_to_use_region` defaults to 10 for WGS:** a region with fewer reads is skipped, its row in
`SAMPLES_QUANTIFICATION_SUMMARY.txt` is `NA`, and the run still exits 0. Check for `NA` rows, and lower the
threshold when the BAM is a small slice. Checked on CRISPResso2 2.3.4 with the CRISPResso2 repo's `tests/`
`smallGenome.fa` (chr9/chr11 slices) and `Both.Cas9.fastq.smallGenome.bam`: FANCF 23 reads, 26.09% Modified;
HEK3 2 reads, `NA` at the default and 50% Modified with `--min_reads_to_use_region 1`.

**Use case:** Validate empirically that an in vivo / clinical-grade edit has minimal off-target activity (combine with GUIDE-seq or CIRCLE-seq predicted sites; randomly chosen regions yield no useful comparison).
