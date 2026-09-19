# Batch Downloads Usage Guide

## Overview

Bulk-download records from NCBI E-utilities efficiently. Encodes the strategy decision (direct fetch vs EPost vs history server vs Datasets CLI), precise rate-limit math, WebEnv lifecycle for long jobs, EPost 200-ID per-call limit, retry/resume design, integrity verification, and the defection rule to NCBI Datasets v2 CLI for genome/gene bulk work.

## Prerequisites

See SKILL.md's Required Setup section for install commands and the `Entrez.email`/`api_key`/`tool` setup code.

## Quick Start

- "Download all human RefSeq mRNAs to a single FASTA, with disk checkpointing so the job can resume"
- "Pull EFetch records for these 5,000 protein accessions; use EPost first to avoid URL-length errors"
- "Download every PubMed abstract for 'CRISPR AND 2024[PDAT]' in MEDLINE format"
- "Compare cost of pulling 100K nucleotide records via EFetch vs NCBI Datasets CLI"
- "Verify the downloaded FASTA has the expected record count and no truncation"

## Example Prompts

### Resumable bulk download

> "Download all RefSeq mRNAs for Homo sapiens to refseq_mrna.fasta. Use the history server. Checkpoint progress to ckpt.json so if my SSH session drops we resume from where we left off, not from zero."

### Large known-ID list

> "I have 5,000 protein accessions in a file. EPost them in chunks of 200 (since EPost's per-call limit is 200), then EFetch by WebEnv/QueryKey in batches of 500."

### Defecting to Datasets CLI

> "I need every RefSeq bacterial genome assembly. Don't loop EFetch -- use 'datasets download genome taxon Bacteria --refseq' from the NCBI Datasets v2 CLI. Show me the equivalent E-utils pipeline so I can see why Datasets is the right tool."

### Production retry / session expiry

> "Build a download with exponential-backoff retry on 429, detection of HTTP-200-with-ERROR-body (WebEnv expired), and automatic re-ESearch + resume from the disk checkpoint."

### Post-download integrity

> "After the download, parse the output FASTA with SeqIO and assert the record count matches what ESearch returned. If not, surface a warning."

## Related Skills

- entrez-search - Build queries before batch-fetching
- entrez-fetch - Single-record fetches
- entrez-link - Chain ELink with neighbor_history for bulk cross-db work
- ncbi-datasets-cli - Modern bulk endpoint for genome and gene data
- sra-data - SRA toolkit for raw sequencing read downloads
- geo-data - GEO supplementary file downloads
