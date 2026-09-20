# BLAST Searches Usage Guide

## Overview

Run remote BLAST searches against NCBI's servers using `Bio.Blast.NCBIWWW`. Encodes the program-vs-database decision matrix, Karlin-Altschul E-value math and why bit-scores beat E-values across databases, the `max_target_seqs` misinterpretation that has biased many published workflows (Shah 2019), composition-based-statistics modes, word-size choice for cross-species vs intra-species DNA, and reproducibility traps with `nt`/`nr`.

## Prerequisites

```bash
pip install biopython
```

No API key required for remote BLAST; one search at a time is the polite cap.

## Quick Start

- "Identify this unknown DNA sequence -- BLASTN against refseq_select for stability, not nt"
- "Find mammalian homologs of this protein in Swiss-Prot using blastp with composition-based statistics mode 2"
- "Run blastp for this short peptide (12 aa) using PAM30 matrix and word size 2"
- "Use dc-megablast for this cross-species mRNA query -- megablast won't find diverged homologs"
- "Search a sequence with hitlist_size=500 and post-filter top 10 by bit-score, avoiding the max_target_seqs trap"

## Example Prompts

### Picking the right program

> "I have an unknown DNA sequence from an environmental sample. Run blastn against refseq_select_rna with word_size=11 and E-value cutoff 1e-10. Don't use megablast -- it requires 28-nt exact match seeds and will miss divergent homologs."

### Reproducibility-safe database choice

> "I'm BLASTing for a publication. Use refseq_select instead of nr/nt because those change daily and aren't reproducible. If I need the broader search, record today's date and archive a frozen copy of the database snapshot."

### Avoiding the max_target_seqs trap

> "Set hitlist_size=500 not 10. The Bio.Blast hitlist_size parameter maps to BLAST+'s max_target_seqs, which is an early-termination heuristic, not 'top N by E-value' (Shah et al. 2019). Then post-filter to the top 10 by bit-score in Python."

### Short-peptide search

> "BLAST this 12-aa proteomics peptide against swissprot. Use matrix='PAM30', word_size=2, expect=1000, composition_based_statistics=3. Default BLOSUM62/word=3 will miss everything for queries this short."

### Cross-database analysis caveat

> "I have results from a blastp vs nr and a blastp vs swissprot for the same query. Don't compare E-values across databases -- they scale with database size. Sort by bit-score for cross-DB comparison."

## What the Agent Will Do

Program choice, database choice, word-size/matrix/CBS selection, the `max_target_seqs`/
`hitlist_size` pattern, `entrez_query` pre-filtering, XML parsing, and bit-score vs. E-value
sorting are all in `SKILL.md`'s Program decision, Database decision, and Code patterns sections --
followed as documented there, including the corrected `gapcosts`/megablast patterns.

## Related Skills

- local-blast - Local BLAST+ for >50 sequences or custom databases
- remote-homology - PSI-BLAST, jackhmmer, HHblits, MMseqs2, DIAMOND, Foldseek for distant homology
- ortholog-inference - RBH, OrthoFinder, OMA for orthology calls
- sequence-io/read-sequences - Load query FASTA files
- entrez-fetch - Fetch full GenBank records for top hits
