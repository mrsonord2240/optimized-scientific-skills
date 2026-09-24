# Entrez Link Usage Guide

## Overview

Navigate between NCBI databases using `Bio.Entrez.elink()`. The skill encodes the most consequential decision in ELink work -- which `linkname` to use for each (`dbfrom`, `db`) pair -- and the auxiliary decisions: which `cmd` variant (`neighbor`, `neighbor_history`, `acheck`, `neighbor_score`), how to handle batches >200 IDs via EPost + history, and how to guard against the asymmetry of link tables.

See SKILL.md's "Required Setup" for install and `Entrez.email`/`api_key` setup, and "Workflow" for the step-by-step agent flow.

## Quick Start

- "Find RefSeq proteins for human BRCA1 (Gene UID 672) using the curated link, not the noisy all-protein link"
- "Get GeneRIF-curated genes mentioned in a given PubMed paper (pubmed_gene_rif)"
- "Discover all available link tables for a gene record"
- "Link 5,000 PubMed UIDs to genes -- the batch is too big for a comma-joined URL"
- "Find SRA runs for BioProject PRJNA123456"

## Example Prompts

### Picking the right linkname

> "For gene UID 672, get the linked proteins using linkname='gene_protein_refseq' so we get the RefSeq isoforms only instead of the larger all-proteins variant."

### Enumerating link options

> "Before I build this pipeline, show me every available linkname for (dbfrom=gene, source-id=672) using cmd='acheck'. I want to see what NCBI exposes and what each link's curation level is."

### Asymmetric round-trip awareness

> "For TP53 (Gene UID 7157), get the curated GeneRIF PubMed citations (gene_pubmed_rif) and compare the count to all PubMed citations (gene_pubmed) -- confirm the curated set is a proper subset, and warn me whenever a curated-vs-all comparison like this could hide the true asymmetry of a round trip through a differently-curated linkname pair."

### Large batch via history server

> "Link 5,000 gene UIDs to RefSeq proteins. The id list is too long for a comma-joined URL, so EPost the UIDs in chunks of 200, then ELink with cmd='neighbor_history' and hand the resulting WebEnv to downstream EFetch."

### BioProject to SRA

> "I have BioProject PRJNA123456. Resolve to a UID, then ELink to SRA to get the run UIDs. Hand off to the sra-data skill to download FASTQ."

See SKILL.md's "Workflow" section for the step-by-step agent flow, "Failure modes" for pitfalls (wrong linkname, empty LinkSetDb, asymmetric round-trips, URL length, per-input indexing), and "Per-database link catalog" for curated-vs-umbrella linkname choices per database pair.

## Related Skills

- entrez-search - Resolve accessions/symbols to UIDs before linking
- entrez-fetch - Retrieve content for linked UIDs
- batch-downloads - Pull large history-server linksets efficiently
- geo-data - GEO-specific links (gds <-> sra, pubmed, bioproject)
- ncbi-datasets-cli - Modern CLI for many gene/genome cross-reference queries
- local-blast / remote-homology - for actual sequence/structure similarity rather than curated ELink relationships
