---
name: bio-alignment-pairwise
description: Perform pairwise sequence alignment using Biopython Bio.Align.PairwiseAligner (Needleman-Wunsch global, Smith-Waterman local, semiglobal). Use when comparing two sequences, finding optimal alignments, scoring similarity, computing percent identity, checking the reverse-complement strand, reproducing EMBOSS needle/water or BLAST scores, and identifying local or global matches between DNA, RNA, or protein sequences.
tool_type: python
primary_tool: Bio.Align
license: MIT
author: GPTomics
---

## Version Compatibility

Checked on Biopython 1.88 (`pip install biopython`), with scores cross-checked against EMBOSS 6.6.0 `needle`/`water`, BLAST+ 2.17.0, parasail 1.3.4, edlib 1.3.9, pywfa 0.5.1, mappy 2.31 and R pwalign 1.2.0 (Bioconductor 3.20).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Pairwise Sequence Alignment

**"Align two sequences"** -> Compute an optimal alignment between a pair of sequences using dynamic programming.
- Python: `PairwiseAligner()` (BioPython Bio.Align)
- CLI: `needle` (global) or `water` (local) from EMBOSS
- R: `pwalign::pairwiseAlignment()` (Bioconductor 3.20 moved it out of Biostrings; `Biostrings::pairwiseAlignment()` still works but warns). Gap parameters use the BLAST convention: see Gap Penalties

Align two sequences using dynamic programming algorithms (Needleman-Wunsch for global, Smith-Waterman for local).

## Required Import

**Goal:** Load modules needed for pairwise alignment operations.

**Approach:** Import the PairwiseAligner class along with sequence and I/O utilities from Biopython.

```python
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio import SeqIO
```

## Core Concepts

| Mode | Algorithm | Use Case |
|------|-----------|----------|
| `global` | Needleman-Wunsch | Full-length alignment, similar-length sequences |
| `local` | Smith-Waterman | Find best matching regions, different-length sequences |
| `global` + free end gaps | Semi-global | Overlap detection, fragment-to-reference alignment |

### Choosing the Right Mode

- **Global**: Both sequences are expected to be homologous over their full length (e.g., two orthologs of similar size). Forces end-to-end alignment.
- **Local**: Conserved domains or motifs within otherwise dissimilar sequences. BLAST uses local alignment internally. Preferred when protein termini are highly divergent (termini accumulate mutations faster than core regions).
- **Semi-global**: One sequence is a fragment or subsequence of the other (e.g., primer to template, read to reference, detecting overlap between shotgun reads). Free end gaps prevent penalizing unaligned flanking regions.

**Common mistake**: Global alignment of sequences with very different lengths forces biologically meaningless terminal gaps. If sequences differ substantially in length, use local or semi-global instead. For pairs beyond a few thousand residues or high-throughput screens, see `references/library-selection.md`.

### DNA vs Protein Alignment

| Scenario | Align As | Rationale |
|----------|----------|-----------|
| Nucleotide identity >70% | DNA | Sufficient signal at nucleotide level |
| Nucleotide identity <70% | Protein | Codon degeneracy masks signal at DNA level; protein alignment is ~3x more sensitive |
| Noncoding sequences (UTRs, intergenic) | DNA | No protein translation possible |
| Coding sequences for dN/dS analysis | Protein first, then back-translate codons (PAL2NAL) | Preserves reading frame for selection analysis |

When in doubt, align at the protein level. It captures functional constraint better because 20 amino acids provide richer signal than 4 nucleotides.

Back-translating a protein alignment onto its CDS (PAL2NAL v14, MAFFT 7.526):

```bash
mafft --auto prot.fa > prot_aln.fa                                # any protein aligner; record IDs must match nuc.fa
pal2nal.pl prot_aln.fa nuc.fa -output fasta > codon_aln.fa        # 3 nucleotide columns per protein column
```

`pal2nal.pl` prints `#--- ERROR: inconsistency between the following pep and nuc seqs ---#` and **exits 0 with an empty output file** when a CDS does not translate to its protein (internal stop, wrong frame). Check the output size, not the exit code: human vs cow HBB CDS gave 441 codon columns (3 x 147); rabbit NM_001314043.1 gave 0 bytes.

## Creating an Aligner

**Goal:** Configure a PairwiseAligner with appropriate scoring for the sequence type.

**Approach:** Instantiate PairwiseAligner with mode, scoring parameters, or a substitution matrix depending on DNA vs protein input.

```python
# Basic aligner with defaults
aligner = PairwiseAligner()

# Configure mode and scoring
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)

# For protein alignment with substitution matrix
from Bio.Align import substitution_matrices
aligner = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'))
```

## Performing Alignments

**"Align two sequences"** -> Compute optimal alignment(s) between a pair of sequences, returning alignment objects or a score.

**Goal:** Align two sequences and retrieve the optimal alignment(s) or score.

**Approach:** Call `aligner.align()` for full alignment objects or `aligner.score()` for score-only (faster for large sequences).

```python
seq1 = Seq('ACCGGTAACGTAG')
seq2 = Seq('ACCGTTAACGAAG')
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)

# Optimal alignments (lazy iterable; see Iterating Over Multiple Alignments)
alignments = aligner.align(seq1, seq2)
print(alignments[0])  # Print first alignment

# Get score only (faster for large sequences)
score = aligner.score(seq1, seq2)
```

## Alignment Output Format

```
target            0 ACCGGTAACGTAG 13
                  0 ||||.|||||.|| 13
query             0 ACCGTTAACGAAG 13
```

## Accessing Alignment Data

**Goal:** Extract alignment properties including score, shape, aligned sequences, and coordinate mappings.

**Approach:** Access alignment object attributes and indexing to retrieve per-sequence aligned strings and coordinate arrays.

```python
alignment = alignments[0]

# Basic properties
print(alignment.score)                    # Alignment score
print(alignment.shape)                    # (num_seqs, alignment_length)
print(len(alignment))                     # Alignment length

# Get aligned sequences with gaps
target_aligned = alignment[0, :]          # First sequence (target) with gaps
query_aligned = alignment[1, :]           # Second sequence (query) with gaps

# Get coordinate mapping
print(alignment.aligned)                  # Array of aligned segment coordinates
print(alignment.coordinates)              # Full coordinate array
```

## Alignment Counts (Identities, Mismatches, Gaps)

**Goal:** Quantify identities, mismatches, and gaps in an alignment to calculate percent identity.

**Approach:** Use the `.counts()` method on the alignment object and derive percent identity from identity and mismatch totals. Percent-identity definitions differ (PID1-4): `references/percent-identity.md`.

```python
alignment = alignments[0]
counts = alignment.counts()

print(f'Identities: {counts.identities}')
print(f'Mismatches: {counts.mismatches}')
print(f'Gaps: {counts.gaps}')

# Calculate percent identity
total_aligned = counts.identities + counts.mismatches
percent_identity = counts.identities / total_aligned * 100
print(f'Percent identity: {percent_identity:.1f}%')
```

## Common Scoring Configurations

### Gap Penalties: Set Them Explicitly, and Mind the Convention

`PairwiseAligner()` with no arguments (Biopython 1.88) uses match_score=1, mismatch_score=0 and open_gap_score=extend_gap_score=-1 (`print(aligner)` shows the resolved set). These defaults are tuned for nothing: combined with BLOSUM62, gaps are nearly free and alignments fragment (HBA vs HBB: 55 gap positions in 37 aligned blocks with the defaults, versus 9 gap positions in 5 blocks at open -11 / extend -1). **Always specify gap penalties explicitly.**

Tools disagree on what "open" means, so the same numbers give different scores:

| Convention | Cost of a gap of length k | Used by |
|---|---|---|
| open + (k-1) * extend | open covers the first gap position | Biopython `open_gap_score`, EMBOSS `needle`/`water`, parasail |
| open + k * extend | open is charged on top of the extension | BLAST+ `-gapopen`, R pwalign `gapOpening` |

Conversion: Biopython `open_gap_score = -(BLAST gapopen + gapextend)`, `extend_gap_score = -gapextend`. So **BLASTP defaults (BLOSUM62, 11/1) are `open_gap_score=-12, extend_gap_score=-1`**, and `-11/-1` is EMBOSS 11/1 (pwalign `gapOpening=10, gapExtension=1`). EMBOSS defaults 10/0.5 are `-10/-0.5`.

Verified on HBA_HUMAN vs HBB_HUMAN, BLOSUM62, local: Biopython -12/-1 = 285 = `blastp -comp_based_stats 0` raw score (same HSP, query 3-141) = pwalign `gapOpening=11, gapExtension=1`; Biopython -11/-1 = 288 = EMBOSS `water -gapopen 11 -gapextend 1` = pwalign `gapOpening=10, gapExtension=1`. Global: -11/-1 = 286 = `needle` 11/1. BLASTP's own default composition-based statistics change the reported score (286 here); use `-comp_based_stats 0` to compare raw scores.

EMBOSS command lines for the same numbers (6.6.0, one sequence per FASTA file; `-auto` suppresses prompts; the matrix file is `EBLOSUM62`, not `BLOSUM62`; DNA uses the default `EDNAFULL`, gaps 10/0.5):

```bash
needle -auto -asequence a.fa -bsequence b.fa -gapopen 11 -gapextend 1 -datafile EBLOSUM62 -outfile needle.txt   # global; "# Score:" line
water  -auto -asequence a.fa -bsequence b.fa -gapopen 11 -gapextend 1 -datafile EBLOSUM62 -outfile water.txt    # local
```

### DNA/RNA Alignment
```python
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
```

### Protein Alignment
```python
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load('BLOSUM62')
aligner = PairwiseAligner(mode='global', substitution_matrix=blosum62, open_gap_score=-11, extend_gap_score=-1)  # EMBOSS 11/1; BLASTP-equivalent is -12/-1; other matrices: references/substitution-matrices.md
```

### Local Alignment (Find Best Region)
```python
aligner = PairwiseAligner(mode='local', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
```

### Semiglobal (Overlap/Fragment Alignment)
```python
# Free end gaps on BOTH sequences (order-independent) -- overlap detection between two reads,
# or a fragment inside a reference
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
aligner.end_gap_score = 0.0
alignment = aligner.align(reference, fragment)[0]   # 20-nt fragment in a 620-nt reference: score 40, target span [[300, 320]]

# Free end gaps on the query (second argument) only. Argument ORDER matters: the fragment must be
# the second argument (align(fragment, reference) scored -279 instead of 40).
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
aligner.open_left_deletion_score = aligner.extend_left_deletion_score = 0
aligner.open_right_deletion_score = aligner.extend_right_deletion_score = 0
alignment = aligner.align(reference, fragment)[0]
```
Older Biopython named these `query_left_open_gap_score`, `query_left_extend_gap_score`, `query_right_open_gap_score`, `query_right_extend_gap_score` (still accepted on 1.88, with a DeprecationWarning).

## Working with SeqRecord Objects

**Goal:** Align sequences loaded from FASTA files rather than hardcoded strings.

**Approach:** Parse SeqRecord objects from a FASTA file and pass their `.seq` attributes to the aligner.

```python
from Bio import SeqIO

records = list(SeqIO.parse('sequences.fasta', 'fasta'))
seq1, seq2 = records[0].seq, records[1].seq

aligner = PairwiseAligner(mode='global', match_score=1, mismatch_score=-1)
alignments = aligner.align(seq1, seq2)
```

## Iterating Over Multiple Alignments

`aligner.align()` returns a lazy iterable. `PairwiseAligner` has no `max_alignments` attribute on Biopython 1.88 (assigning it raises AttributeError). `len(alignments)` counts every optimal alignment and raises `OverflowError` when there are too many (e.g. zero gap penalties on repetitive input), so cap with `itertools.islice` when the count is not needed:

```python
from itertools import islice

for i, alignment in enumerate(islice(alignments, 5)):
    print(f'Alignment {i+1}: score={alignment.score}')
```

## Quick Reference: Scoring Parameters

| Parameter | Description | Typical DNA | Typical Protein |
|-----------|-------------|-------------|-----------------|
| `match_score` | Score for identical bases | 1-2 | Use matrix |
| `mismatch_score` | Penalty for mismatches | -1 to -3 | Use matrix |
| `open_gap_score` | Cost to start a gap | -5 to -15 | -10 to -12 |
| `extend_gap_score` | Cost per gap extension | -0.5 to -2 | -0.5 to -1 |
| `substitution_matrix` | Scoring matrix | N/A | BLOSUM62 |

## Input Checks

Verified on Biopython 1.88 with BLOSUM62 / NUC.4.4 aligners:
- **Case and whitespace**: lowercase, a trailing newline, or `J`/`U` residues raise `ValueError` with a substitution matrix (NUC.4.4 also rejects `U`). Use `str(seq).strip().upper()`. Match/mismatch aligners do not raise: `'acgt'` vs `'ACGT'` scores as four mismatches, and `U` vs `T` as a mismatch (convert RNA to DNA first).
- **Empty sequences** raise `ValueError: sequence has zero length`. Check lengths before aligning.
- **Accepted silently**: `*` (stop), `X`/`B`/`Z`, and a `SeqRecord` (scores identically to its `.seq`; pass `.seq` explicitly).
- **Strand**: alignment is strand-specific. A reverse-complemented query scores far lower (a random 30-nt exact match scored 60 in local mode; its reverse complement scored a median 18, never above 29, over 300 random pairs; palindromic or low-complexity queries score higher). For DNA of unknown orientation, score both `seq` and `seq.reverse_complement()` and keep the higher.
- **Coding sequences**: check for internal stops (`'*' in str(seq.translate()).rstrip('*')`) before aligning or back-translating; a CDS with an internal stop (e.g. RefSeq NM_001314043.1, rabbit HBB2) breaks codon-aware tools such as PAL2NAL.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `OverflowError` | `len(alignments)` with too many optimal alignments | Iterate with `itertools.islice(alignments, n)` |
| Low scores | Wrong scoring scheme | Use substitution matrix for proteins |
| No alignments in local mode | Scores all negative | Ensure `match_score` > 0 |
| `ValueError: sequence contains letters not in the alphabet` / `zero length` | See Input Checks | Clean the input first |

## Reference Files

| File | Read when |
|------|-----------|
| `references/library-selection.md` | pairs beyond a few thousand residues, throughput, parasail / edlib / pywfa / mappy snippets, saturation |
| `references/substitution-matrices.md` | choosing BLOSUM/PAM/NUC by divergence, affine gap rationale |
| `references/percent-identity.md` | reporting percent identity (PID1-4 definitions) |
| `references/significance.md` | bit score / E-value, empirical p-values, compositional bias |
| `references/when-not-appropriate.md` | low identity (<40%), database-scale search, MMseqs2 / jackhmmer / HHsearch command lines, repeats |
| `references/alignment-export.md` | `.substitutions` counts, FASTA/Clustal/PSL/SAM export |

## Related Skills

- alignment/multiple-alignment - Align three or more sequences with MAFFT, MUSCLE5, ClustalOmega
- alignment/alignment-io - Save alignments to files in various formats
- alignment/msa-parsing - Work with multiple sequence alignments
- alignment/msa-statistics - Calculate identity, similarity metrics
- alignment/structural-alignment - Twilight-zone alternative when sequence signal fails (Foldseek, TM-align, pLM aligners)
- alignment/alignment-trimming - Remove unreliable columns post-alignment
- sequence-manipulation/motif-search - Pattern matching in sequences
