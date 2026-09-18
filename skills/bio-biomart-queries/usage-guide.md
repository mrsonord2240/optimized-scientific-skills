# BioMart Queries Usage Guide

## Overview

Bulk-query Ensembl BioMart for cross-database ID mapping, coordinate tables, and ortholog wide tables -- for anything expected to return >5,000 rows, where looping Ensembl REST would be slow. See SKILL.md's "Decision matrix" for the BioMart-vs-REST call, and "Installation" for setup.

## Quick Start

- "Convert 5,000 Ensembl Gene IDs to HGNC + RefSeq + UniProt in one query"
- "Pull all protein-coding genes on chr17 with coordinates and biotype"
- "Get a wide ortholog table: human Ensembl Gene ID, mouse ortholog, zebrafish ortholog"
- "Fetch GO term annotations for a list of genes (long format)"
- "Pin BioMart to release 110 for reproducibility"

## Example Prompts

### Bulk ID mapping

> "I have 8,000 Ensembl Gene IDs. Convert them to HGNC symbol, NCBI Entrez Gene ID, RefSeq mRNA accessions, and Swiss-Prot UniProt accessions. Use one pybiomart query against hsapiens_gene_ensembl with filters={'ensembl_gene_id': [...]}. Don't loop Ensembl REST -- that's 8,000 sequential calls."

### Coordinate table

> "Pull all protein-coding genes on chromosome 17 with: ensembl_gene_id, external_gene_name, start, end, strand, biotype. One BioMart query."

### Ortholog wide table

> "Build a wide table: for every human protein-coding gene, the mouse ortholog Ensembl ID and the zebrafish ortholog Ensembl ID, both filtered to ortholog_one2one only. Use mmusculus_homolog_ensembl_gene + drerio_homolog_ensembl_gene attributes."

### GO annotations long format

> "Get GO term annotations for [TP53, BRCA1, MYC, EGFR]. Returns long format -- one row per (gene, GO term)."

### Version pinning

> "I'm writing a paper. Use R biomaRt with useEnsembl(version=110) so the query reproduces in 2030. Don't use the default which follows the current release."

## Related Skills

- ensembl-rest - Per-record Ensembl queries (complement to BioMart)
- ortholog-inference - Compara orthologs with confidence semantics
- uniprot-access - Alternative ID-mapping via UniProt
- ncbi-datasets-cli - NCBI-side bulk path
- entrez-search - NCBI alternative for non-Ensembl queries
