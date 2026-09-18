# Local BLAST Usage Guide

## Overview

Build local BLAST+ databases and run searches using the BLAST+ CLI. Encodes BLAST v5 vs v4 database format and taxonomy filtering, `-task` taxonomy (megablast/dc-megablast/blastn/blastn-short and the blastp variants), soft vs hard masking, thread saturation past ~16 threads, output-field semantics (pident vs qcovs vs qcovhsp), the `max_target_seqs` misuse (Shah 2019), the `-parse_seqids` requirement for downstream extraction, and the size/reproducibility trap of `nt`/`nr`.

## Prerequisites

```bash
conda install -c bioconda blast
blastn -version       # 2.15+ expected
update_blastdb.pl --showall pretty | head
```

For large workflows (>100,000 queries), consider DIAMOND or MMseqs2 instead -- 100-10,000x faster at comparable sensitivity. See `remote-homology`.

## Quick Start

- "Build a v5 BLAST database from my reference FASTA with -parse_seqids so I can extract hits later"
- "Run blastp against my custom database with 8 threads and hitlist_size=500"
- "Download refseq_select_rna with update_blastdb.pl (5 GB) instead of nt (250 GB)"
- "Use dc-megablast for cross-species cDNA; default megablast misses divergent hits"
- "Filter blastn results to top hit per query by bit-score, with qcovs >= 0.8"

## Example Prompts

### Building a custom DB correctly

> "Build a protein BLAST database from reference_proteins.fasta. Use -blastdb_version 5 and -parse_seqids so blastdbcmd can extract sequences by accession later, and -hash_index for faster lookups."

### Picking the right -task

> "I'm doing cross-species DNA homology between mouse and human. Use -task dc-megablast (discontiguous, sensitive). Don't use default megablast -- that needs 28-nt exact-match seeds and misses cross-species hits."

### Avoiding the size shock

> "Download NCBI's curated set for cross-species homology. Use refseq_select_rna (~5 GB) via update_blastdb.pl, not nt (~250 GB). nt also isn't reproducible -- it changes daily."

### Thread saturation awareness

> "I have a 32-core box. Don't set -num_threads 32 -- BLAST is I/O bound past ~16 threads. Either cap at 16 threads, or split the query FASTA into 4 chunks and run 4 parallel blastn processes with -num_threads 8 each."

### Avoiding max_target_seqs misuse

> "Set -max_target_seqs 500 and post-filter to top 10 by bit-score. Don't set -max_target_seqs 10 -- it's an early-termination heuristic that biases what's returned (Shah et al. 2019 Bioinformatics 35:1613)."

### Hit extraction

> "After blastp finishes, use blastdbcmd -db ref_prot_db -entry_batch hit_accessions.txt to pull the subject FASTA sequences for downstream analysis."

## What the Agent Will Do

1. Verify BLAST+ install and version with `blastn -version`.
2. Build databases with v5 format + `-parse_seqids` + `-hash_index`.
3. Pick `-task` deliberately: `megablast` for >=95% identity, `dc-megablast` for cross-species, `blastn` for general, `blastn-short` for <50 nt.
4. Use `refseq_select` over `nt`/`nr` for reproducibility unless explicitly needed otherwise.
5. Set `-num_threads` between 8 and 16; split input for higher parallelism.
6. Use `-outfmt "6 qseqid sseqid pident length qcovs qcovhsp evalue bitscore staxids sscinames stitle"` for analyzable tabular.
7. Set `-max_target_seqs 500+` and post-filter; cite Shah 2019 if the user asks why.
8. Use `-taxids` / `-taxidlist` for taxonomy filtering (v5 only); fall back to entrez_query for v4 (remote).
9. Extract hit sequences with `blastdbcmd` only if DB was built with `-parse_seqids`.

## Tips

See SKILL.md's "Database format: v5 vs v4" (taxonomy filtering prerequisites), "Soft vs hard masking",
"Output format reference" (`pident` vs `qcovs`), "Thread scaling" (DIAMOND/MMseqs2 for >100K queries),
"Code patterns" (RBH paralog caveat), and "Practice boundaries" for the details an agent needs --
kept in one place there instead of restated here.

## Related Skills

- blast-searches - Remote BLAST against NCBI servers
- remote-homology - PSI-BLAST, jackhmmer, HHblits, MMseqs2, DIAMOND, Foldseek
- ortholog-inference - Principled orthology via RBH/OrthoFinder/OMA
- sequence-io/read-sequences - Load FASTA inputs/outputs
- batch-downloads - Download reference FASTA before makeblastdb
