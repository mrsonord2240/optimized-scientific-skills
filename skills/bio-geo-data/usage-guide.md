# GEO Data Usage Guide

## Overview

Query NCBI GEO (and EMBL-EBI's BioStudies/ArrayExpress mirror) for expression datasets. Encodes the SuperSeries trap (a GSE may wrap multiple SubSeries on different platforms; default download mixes them), the series-matrix normalization-trust caveat (submitter-normalized; re-process from raw for reproducibility), processed-vs-raw decision (Affymetrix CEL + locally-run RMA; RNA-seq via SRA FASTQ), GEOparse vs GEOquery trade-off, GEOmetadb deprecation (2020), and ArrayExpress migration to BioStudies (2020).

## Prerequisites

See SKILL.md's "Required Setup" section.

## Quick Start

- "Search GEO for human single-cell RNA-seq from 2024; flag any SuperSeries"
- "Detect whether GSE122288 is a SuperSeries before downloading; if so, list its SubSeries"
- "Resolve GSE123456 to SRA run accessions via pysradb; hand off to sra-data for FASTQ download"
- "Download the series matrix for GSE123456 and tell me what 'normalization' the submitter applied"
- "Pull supplementary CEL files for GSE12345 via R's GEOquery -- GEOparse has been flakey on suppl downloads"

## Example Prompts

### Detecting SuperSeries before pulling

> "Before I download GSE122288 as a single experiment, check its SOFT family file for !Series_relation. If it's a SuperSeries, list each SubSeries and recommend processing them independently to avoid mixing platforms."

### Processed-vs-raw decision

> "I want expression values for GSE123456. Check whether the technology is Affymetrix (-> download CEL files and re-do RMA) or RNA-seq (-> resolve to SRA and re-quantify with Salmon). Don't trust the submitter's series-matrix values for downstream stats."

### GEO -> SRA linkage

> "Find SRA accessions for GSE123456 using pysradb (more reliable than the gds -> sra ELink). Return a list of SRR run IDs to hand off to the sra-data skill."

### Cross-reference from publication

> "Find all GEO datasets cited in PMID 35412348 via pubmed -> gds ELink. Summarize each with title, sample count, platform, and SuperSeries status."

### Submitter-data-processing audit

> "Download the series matrix for GSE123456 and dump every unique value of !Sample_data_processing. Tell me whether the matrix is raw counts, log-CPM, VST, or RMA-normalized -- the answer determines whether I can use it as-is."

## What the Agent Will Do

See SKILL.md's "Workflow" section.

## Related Skills

- entrez-search - General gds search
- entrez-link - gds <-> sra / pubmed / bioproject links
- sra-data - Pull raw FASTQ from GEO-linked SRA
- expression-matrix/normalization - Re-normalize raw expression data
- rna-quantification/alignment-free-quant - Salmon/kallisto re-quantification
- ensembl-rest - Cross-reference Ensembl IDs in series matrices
