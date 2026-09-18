# SRA Data Usage Guide

## Overview

Download raw sequencing reads from NCBI SRA. Encodes the source-of-truth decision (SRA-direct vs ENA mirror vs AWS/GCP STRIDES cloud), prefetch `--max-size` trap (silent 20 GB skip), fasterq-dump uncompressed-scratch trap (~3x final size), `--include-technical` for 10x single-cell records, MD5 validation, accession hierarchy navigation (SRR/SRX/SRS/SRP/PRJNA), pysradb for metadata, the controlled-access (dbGaP) boundary, and Aspera deprecation post-2019.

## Prerequisites

See SKILL.md's "Required Setup" section for sra-tools/pysradb install commands, vdb-config cache setup, and the AWS CLI/STRIDES note.

## Quick Start

- "Download SRR12345678 as paired-end FASTQ via the ENA mirror with MD5 verification"
- "Resolve GSE123456 to its SRR run accessions using pysradb, then download with prefetch"
- "Download a 10x single-cell record -- include technical reads (barcodes/UMIs/indexes)"
- "Pull SRA data from AWS Open Data (STRIDES) into an EC2 instance in us-east-1"
- "Validate downloaded FASTQ against ENA's md5; fail loudly on mismatch"

## Example Prompts

### Default: ENA mirror

> "Download SRR12345678 via the ENA mirror (https://ftp.sra.ebi.ac.uk/). It's typically faster than SRA-direct because the FASTQ is pre-compressed -- no SRA->FASTQ conversion needed. Get the md5 from the ENA portal API and verify."

### Resolving GSE to SRR

> "I have GSE123456 from a GEO paper. Use pysradb gse_to_srp -> srp_to_srr to get all SRR run accessions. Then batch-download them via ENA."

### Large prefetch without silent skip

> "Use prefetch with --max-size 200G explicitly -- the default 20 GB silently skips larger runs. After prefetch, run vdb-validate and then fasterq-dump."

### 10x single-cell

> "This is a 10x v3 record. fasterq-dump with --include-technical --split-files. Expect 3 files: R1 (28-bp barcode+UMI), R2 (cDNA), I1 (sample index). The default (without --include-technical) only gives R2 which is useless for CellRanger."

### Cloud-native pipeline

> "I'm running a Nextflow pipeline on AWS Batch in us-east-1. Pull SRA data from s3://sra-pub-run-odp/ with --no-sign-request. Same-region transfer is free."

### Scratch disk management

> "fasterq-dump writes uncompressed FASTQ to scratch (~3x final compressed size). My scratch dir has 500 GB free; the run is 200 GB compressed. That's tight -- use fastq-dump --gzip instead, which writes compressed in-place."

For the agent's decision process (source selection, `--max-size` sizing, pigz/gzip fallback, the ENA-vs-SRA-direct read-count caveat, and the controlled-access/dbGaP boundary), see SKILL.md's decision matrix, "prefetch and the `--max-size` trap", "Controlled-access (dbGaP) data" and "Failure modes" sections -- summarized once there, not repeated here.

## Related Skills

- entrez-search - Search SRA db for accessions
- geo-data - GEO Series often link to SRA via ELink
- read-qc/quality-reports - QC the downloaded FASTQ
- read-qc/fastp-workflow - Adapter trim downloaded FASTQ
- read-alignment/bwa-alignment - Align downloaded reads
- ncbi-datasets-cli - Modern bulk path for genomes (NOT for raw reads)
