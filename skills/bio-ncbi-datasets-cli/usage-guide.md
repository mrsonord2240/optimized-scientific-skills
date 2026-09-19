# NCBI Datasets CLI Usage Guide

## Overview

Use the NCBI Datasets v2 CLI (launched 2023, current as of 2024) for genome and gene-centric bulk workflows. Encodes the defection rule (use Datasets for genome/gene; E-utilities for PubMed/SRA/custom queries), the `--dehydrated` flag for cloud-friendly parallel pulls, JSON-lines output + `dataformat` conversion to TSV, automatic MD5 verification, and the choice between `--reference` (one per species) vs full set.

See SKILL.md's Installation and Version Compatibility sections for setup and the installed-version
check (this doc was checked against 18.37.0).

## Quick Start

- "Download human reference genome (GCF_000001405.40) with genome + GFF3 + protein + CDS"
- "Pull every reference bacterial genome from RefSeq; use --dehydrated + aria2c for parallel transfer"
- "Get a TSV of all reference E. coli assemblies with N50 and assembly level"
- "Find NCBI ortholog set for human BRCA1 across mammals (one rep per species)"
- "Download SARS-CoV-2 genome assemblies released after 2024-01-01"

## Example Prompts

### Single-assembly download

> "Download the human GRCh38 reference assembly (GCF_000001405.40) via Datasets CLI. Include genome,gff3,gtf,protein,cds. Datasets verifies MD5 automatically -- no manual checksum step needed."

### Bulk download via --dehydrated

> "Pull every RefSeq reference bacterial genome with annotation. The serial path takes hours; use --dehydrated to get a fetch.txt of URLs, then aria2c --max-concurrent-downloads=8 for parallel transfer."

### Bulk metadata via dataformat tsv

> "Get a TSV of all reference Salmonella enterica assemblies released since 2024-01-01 with: accession, organism-name, assembly level, scaffold N50, release date. Use datasets summary + dataformat tsv genome --fields=... (check `dataformat tsv genome --help` for the current field names -- see SKILL.md)."

### NCBI ortholog set

> "Get NCBI's ortholog set for human BRCA1 across Mammalia. Use datasets summary gene symbol BRCA1 --ortholog Mammalia (see SKILL.md's 'Gene metadata across species' section -- --taxon alone can't do a cross-species query, --ortholog is the mechanism). Note: this is the simple single-representative-per-species view; for tree-reconciled orthology with co-orthologs use Ensembl Compara via ortholog-inference."

### Datasets vs E-utilities decision

> "I need 100 reference genomes. Don't loop EFetch -- use datasets download genome accession ... It's 5-50x faster, handles checksums, and parallelizes within one ZIP. For PubMed or SRA reads, stay with E-utilities -- Datasets doesn't cover those."

The agent decision flow (scope check, `summary` vs `download`, `dataformat` field lookup,
`--dehydrated` threshold, `--reference`, `--api-key`) and all tips/gotchas live in SKILL.md -- see
"What's in scope", "Subcommand taxonomy", "Key parameters", "When to use --dehydrated", and
"Failure modes" / "Common errors" there. This file stays limited to prompts and related Skills so
there's one place to fix a command, not two.

## Related Skills

- entrez-search - PubMed and non-genome queries
- entrez-fetch - Single-record fetches outside genome/gene scope
- batch-downloads - Bulk E-utilities for non-genome data
- sra-data - Raw sequencing reads (NOT covered by Datasets)
- ensembl-rest - Ensembl REST as alternative
- ortholog-inference - Compara/OMA/OrthoDB for tree-aware orthology
