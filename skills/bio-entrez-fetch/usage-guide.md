# Entrez Fetch Usage Guide

## Overview

Retrieve full records (`EFetch`) or lightweight document summaries (`ESummary`) from NCBI databases using Biopython's `Bio.Entrez`. The skill encodes the rettype/retmode decision matrix per database, ESummary-vs-EFetch triage for bulk metadata work, GI deprecation (records after 2017 have only accession.version), the `gbwithparts` trap for WGS assemblies, XML schema drift, and accession-versioning for reproducible analyses.

## Prerequisites

```bash
pip install biopython
```

See SKILL.md's "Required Setup" section for the `Entrez.email`/`api_key` setup NCBI requires.

## Quick Start

- "Fetch the GenBank record for NM_007294.4 and parse the CDS coordinates"
- "Get organism + length for these 5,000 nucleotide UIDs without downloading sequence"
- "Download the CDS-translated proteins from the E. coli K-12 reference genome in one EFetch call"
- "Pull the full PubMed XML for PMID 35412348 including MeSH terms"
- "Convert a list of SRA UIDs to SRR run accessions with library size metrics"

## Example Prompts

### Choosing ESummary over EFetch

> "I have 8,000 nucleotide UIDs. I only need the organism, accession.version, and sequence length. Use ESummary in chunks of 500, not EFetch -- it's an order of magnitude cheaper."

### Reproducibility via versioned accessions

> "Fetch GenBank for NM_007294 but lock the result to whatever specific .version is current today. Save the accession.version so reruns next year fetch the same record content."

### Handling WGS records correctly

> "Fetch the genome assembly for WGS accession ABFD01000000 -- use rettype='gbwithparts' so the CONTIG sequences are inlined, otherwise the record comes back with no actual sequence."

### One-shot CDS extraction

> "Download all CDS translations from RefSeq NC_000913.3 (E. coli K-12) using rettype='fasta_cds_aa'. Don't walk the GenBank features manually -- let NCBI do the extraction server-side."

### SRA metadata conversion

> "Convert these SRA UIDs to SRR run accessions plus Bases/Spots/AvgLength metrics using EFetch with rettype='runinfo'. Parse the CSV and return as a DataFrame."

For the agent's decision process, the rettype/retmode decision matrix, GI-deprecation rules, the `gbwithparts` WGS trap, XML schema-drift guardrails, and rate-limit/chunking mechanics, see SKILL.md -- summarized once there, not repeated here.

## Related Skills

- entrez-search - Find UIDs to fetch
- entrez-link - Cross-database navigation via ELink
- batch-downloads - History-server pipelines for large fetches
- ncbi-datasets-cli - Modern CLI for genome and gene metadata
- sequence-io/read-sequences - Parse downloaded FASTA/GenBank locally
