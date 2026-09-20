---
name: bio-alignment-pairwise
description: Perform pairwise sequence alignment using Biopython Bio.Align.PairwiseAligner. Use when comparing two sequences, finding optimal alignments, scoring similarity, and identifying local or global matches between DNA, RNA, or protein sequences.
tool_type: python
primary_tool: Bio.Align
license: MIT
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

## Pairwise Library Selection

`Bio.Align.PairwiseAligner` is the right default for interactive use, scripting, and pair sizes up to a few thousand residues, but it is not the fastest or most scalable option. For high-throughput screens, very long sequences, or production pipelines, switch to a SIMD-accelerated or specialised library.

| Library | Speed vs Bio.Align | Alphabet | Scoring | Vectorization | When to use |
|---------|-------------------|----------|---------|---------------|-------------|
| `Bio.Align.PairwiseAligner` (BioPython) | 1x baseline | DNA / RNA / protein | Matrix + affine | C-backed Gotoh | Default, <10 kb pairs, interactive use |
| `parasail` (Daily 2016 BMC Bioinf) | 2-10x (300 nt), ~3-8x (20 kb); measured | DNA / protein | Matrix + affine | SSE / AVX SIMD | High-throughput SW or NW; benchmark loops |
| `edlib` (Sosic & Sikic 2017 Bioinf) | 15-26x (300 nt), 400x+ (20 kb); measured | any alphabet | Edit distance only | Bit-parallel Myers | Read mapping, k-mer search, primer placement |
| `pywfa` / WFA2 (Marco-Sola 2021 Bioinformatics 37:456; BiWFA: Marco-Sola 2023 Bioinformatics 39:btad074) | Best for low-divergence | DNA | Matrix + affine | Wavefront, O(s) memory | Long, near-identical sequences (>10 kb, <5% diverged) |
| `mappy` / minimap2 (Li 2018 Bioinf) | Production reads-to-genome | DNA | Chain + base-level | k-mer chain | Long-read mapping, splice-aware DNA |
| `Bio.pairwise2` | DEPRECATED | -- | -- | -- | Migrate to `PairwiseAligner` (deprecated in BioPython 1.80; not yet removed; migrate proactively) |
| EMBOSS `needle` / `water` | ~Bio.Align | DNA / protein | Matrix + affine | None | Reproducibility, audit trails (fixed, documented default parameters) |

Speed numbers for parasail and edlib were measured (Windows, score only, 5% divergence for 300 nt pairs, 3% for 20 kb; timings vary run to run); the rest are literature figures. Benchmark on representative inputs before committing. Critical caveats:
- **WFA / BiWFA**: 10-100x faster than Gotoh below 5% divergence; above ~10% it converges to Gotoh complexity. Right tool for PacBio HiFi self-similarity or assembly-vs-reference; not for distant homologs.
- **edlib**: measured 15x vs Biopython Levenshtein and 26x vs affine scoring on 1000 pairs of 300 nt, 400x+ on a 20 kb pair; the speedup shrinks as divergence grows (literature: ~64x above ~50% divergence). It returns edit distance only, so for high-divergence DNA (<70% nucleotide identity) prefer parasail's SIMD score-only mode.
- **parasail**: SIMD only realises its advantage on long sequences in amortised batch loops. **Fixed-width variants silently saturate**: `nw_striped_16` on a 20 kb pair returned score 0 with `.saturated == True` (true score 35565). Use the `*_sat` variants (they widen 8 -> 16 -> 32 bit) and check `.saturated`.

Verified snippets (parasail's open/extend use the Biopython convention, so `10, 1` = `open_gap_score=-10, extend_gap_score=-1`):

```python
import parasail, edlib
mat = parasail.matrix_create('ACGT', 2, -1)                  # protein: parasail.blosum62
r = parasail.nw_striped_sat(query, target, 10, 1, mat)       # sw_striped_sat for local
assert not r.saturated
print(r.score)   # == PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-1).score(target, query)

# edit distance == -score of PairwiseAligner(match 0, mismatch -1, all gaps -1); mode 'HW' = query anywhere in target
d = edlib.align(query, target, mode='NW', task='distance')['editDistance']
```

`pywfa` and `mappy` (tested on Linux; their Windows builds failed):

```python
from pywfa import WavefrontAligner
x, o, e = 4, 6, 2                      # WFA2 penalties: a gap of length L costs o + e*L
r = WavefrontAligner(target, span='end-to-end', mismatch=x, gap_opening=o, gap_extension=e)(query)
# r.score == PairwiseAligner(mode='global', match_score=0, mismatch_score=-x, open_gap_score=-(o+e), extend_gap_score=-e).score(target, query)

import mappy
for h in mappy.Aligner('ref.fa', preset='map-ont').map(read):   # presets: sr, map-hifi, asm5, ...
    print(h.ctg, h.r_st, h.r_en, h.mapq, h.cigar_str)
```

When uncertain which algorithm Biopython's aligner selected internally, inspect `aligner.algorithm` after configuration -- it returns the resolved variant ("Needleman-Wunsch", "Smith-Waterman", "Gotoh global alignment algorithm", "Gotoh local alignment algorithm", "Waterman-Smith-Beyer global alignment algorithm", "Waterman-Smith-Beyer local alignment algorithm") for deterministic auditing.

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

**Common mistake**: Global alignment of sequences with very different lengths forces biologically meaningless terminal gaps. If sequences differ substantially in length, use local or semi-global instead.

### DNA vs Protein Alignment

| Scenario | Align As | Rationale |
|----------|----------|-----------|
| Nucleotide identity >70% | DNA | Sufficient signal at nucleotide level |
| Nucleotide identity <70% | Protein | Codon degeneracy masks signal at DNA level; protein alignment is ~3x more sensitive |
| Noncoding sequences (UTRs, intergenic) | DNA | No protein translation possible |
| Coding sequences for dN/dS analysis | Protein first, then back-translate codons (PAL2NAL) | Preserves reading frame for selection analysis |

When in doubt, align at the protein level. It captures functional constraint better because 20 amino acids provide richer signal than 4 nucleotides.

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

**Approach:** Use the `.counts()` method on the alignment object and derive percent identity from identity and mismatch totals.

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

### DNA/RNA Alignment
```python
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
```

### Protein Alignment
```python
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load('BLOSUM62')
aligner = PairwiseAligner(mode='global', substitution_matrix=blosum62, open_gap_score=-11, extend_gap_score=-1)  # EMBOSS 11/1; BLASTP-equivalent is -12/-1
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

## Substitution Matrix Selection

**Goal:** Select the appropriate substitution matrix based on expected sequence divergence.

**Approach:** Match matrix to divergence level. BLOSUM and PAM number in **opposite directions**: higher BLOSUM = closer sequences; higher PAM = more distant sequences.

| Divergence Level | BLOSUM | PAM | When To Use |
|-----------------|--------|-----|-------------|
| Very close (<20% divergence) | BLOSUM80, BLOSUM90 | PAM30 | Recently duplicated genes, strain comparison |
| Moderate | BLOSUM62 (default) | PAM120 | General-purpose, most analyses |
| Distant (>50% divergence) | BLOSUM45, BLOSUM50 | PAM250 | Remote homology detection |

**BLOSUM62 is the universal default** (used by BLAST, most alignment tools). When in doubt, use BLOSUM62. Switch to BLOSUM80 for very similar proteins or BLOSUM45 for distant homologs.

**DNA matrices**: `NUC.4.4` (match=+5, mismatch=-4) accepts IUPAC ambiguity codes but scores partial matches weakly (`R` vs `A` = +1, not +5; check with `print(substitution_matrices.load('NUC.4.4'))` before relying on ambiguity-aware scoring). `HOXD70` is tuned for human-mouse whole-genome alignment from noncoding regions.

```python
from Bio.Align import substitution_matrices
print(substitution_matrices.load())  # List all 30 available matrices

blosum62 = substitution_matrices.load('BLOSUM62')  # General protein (default)
blosum80 = substitution_matrices.load('BLOSUM80')  # Close homologs
blosum45 = substitution_matrices.load('BLOSUM45')  # Distant homologs
nuc44 = substitution_matrices.load('NUC.4.4')      # DNA with IUPAC support
```

### Affine Gap Penalties: Biological Rationale

Gap penalties control how gaps (insertions/deletions) are scored. The **affine model** (`penalty = open + extend * (L-1)`) is almost always preferred over linear because it reflects indel biology: a DNA break introduces the first gap (costly), but extending an existing gap is mechanistically easier (less costly). This models the observation that indels in real sequences tend to occur as single contiguous events.

Typical values with BLOSUM62: `open_gap_score=-12, extend_gap_score=-1` (BLASTP defaults; see Gap Penalties for the convention) or `-11/-1` (EMBOSS). Setting gap open equal to gap extend (linear model) over-penalizes long indels and under-penalizes scattered single-residue gaps, producing biologically unrealistic alignments.

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

## Substitution Matrix from Alignment

**Goal:** Extract observed substitution frequencies from a completed alignment.

**Approach:** Access the `.substitutions` property to get a matrix of observed base/residue substitution counts.

```python
alignment = alignments[0]
substitutions = alignment.substitutions

# View as array (rows=target, cols=query)
print(substitutions)

# Access specific substitution counts
# substitutions['A', 'T'] gives count of A aligned to T
```

## Export Alignment to Different Formats

**Goal:** Convert an alignment to standard bioinformatics file formats for downstream tools.

**Approach:** Use Python's `format()` function with format specifiers (fasta, clustal, psl, sam) on the alignment object.

```python
alignment = alignments[0]

# Various output formats
print(format(alignment, 'fasta'))     # FASTA format
print(format(alignment, 'clustal'))   # Clustal format
print(format(alignment, 'psl'))       # PSL format (BLAT)
print(format(alignment, 'sam'))       # SAM format
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
- **Strand**: alignment is strand-specific. A reverse-complemented query scores near zero (a 30-nt exact match scored 60, reverse complement 0). For DNA of unknown orientation, score both `seq` and `seq.reverse_complement()` and keep the higher.
- **Coding sequences**: check for internal stops (`'*' in str(seq.translate()).rstrip('*')`) before aligning or back-translating; a CDS with an internal stop (e.g. RefSeq NM_001314043.1, rabbit HBB2) breaks codon-aware tools such as PAL2NAL.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `OverflowError` | `len(alignments)` with too many optimal alignments | Iterate with `itertools.islice(alignments, n)` |
| Low scores | Wrong scoring scheme | Use substitution matrix for proteins |
| No alignments in local mode | Scores all negative | Ensure `match_score` > 0 |
| `ValueError: sequence contains letters not in the alphabet` / `zero length` | See Input Checks | Clean the input first |

## Percent Identity: Definitions Matter

There are four common ways to calculate percent identity from the same alignment, producing different values:

| Method | Denominator | Best For |
|--------|-------------|----------|
| PID1 | Aligned positions + internal gaps | Gap-aware, conservative |
| PID2 | Aligned residue pairs (excluding gaps) | Always highest value |
| PID3 | Shorter sequence length | Length-normalized |
| PID4 | Mean sequence length | Best correlation with structural similarity |

**Practical impact**: Up to 11.5% difference between methods on a single alignment. Combined with different alignment algorithms, variation reaches 22%. Always report which method was used. The `counts()` method above uses aligned non-gap positions (similar to PID2).

## Statistical Significance: Karlin-Altschul

Use bit score (database-size-independent) and E-value (expected chance hits) to interpret raw alignment scores; never compare raw scores across scoring schemes. For non-default gap penalties or non-protein/DNA alphabets, generate empirical p-values via sequence shuffling instead of trusting the formula (`examples/empirical_pvalue.py`; for DNA, use the dinucleotide shuffle of Altschul & Erickson 1985 via `ushuffle`).

| Bit score | E-value (typical 1e6 db) | Interpretation |
|-----------|--------------------------|----------------|
| > 200 | < 1e-50 | Essentially certain homology |
| 50-200 | 1e-5 to 1e-50 | Likely homology |
| 30-50 | 0.01 to 1e-5 | Possible homology; verify with profile methods |
| < 30 | > 0.01 | Suspect; not significant |

Karlin-Altschul lambda/K are calibrated empirically per (matrix, gap-open, gap-extend) tuple; switching gap penalties without recalibrating produces wrong E-values.

**Mask compositional bias before scoring.** Low-complexity regions (poly-Q, proline-rich, leucine-rich TM helices) inflate raw scores. Pre-filter with SEG for protein or DUST for DNA, then run BLAST with `-comp_based_stats 2` (modern default; Yu & Altschul 2005 Bioinf, conditional) or `-comp_based_stats 3` for repeat-rich queries (Yu & Altschul 2005 unconditional; Schaffer et al 2001 NAR introduced the original level-1 statistics). SEG pre-filtering is complementary to `-comp_based_stats`, not redundant.

## When Alignment Is NOT Appropriate

Pick the alignment method by protein identity; below 15% identity DP alignments are statistically indistinguishable from random pairings, so escape to structure or pLM tools.

| Identity (protein) | Recommended approach |
|-------------------|---------------------|
| >= 40% | Any DP aligner; Bio.Align or BLAST is sufficient |
| 25-40% | Use sensitive iterative methods (MMseqs2 iterative, BLASTP with composition-based statistics) |
| 15-25% | Profile-profile (HHsearch, HMMER `phmmer`/`jackhmmer`) |
| < 15% | Structural alignment (Foldseek, TM-align, US-align) or pLM embeddings (TM-Vec, ESM-2 + cosine) -- see structural-alignment |

Length amplifies the signal: 30% identity over 200 residues is far more reliable than 30% over 50.

**When pairwise becomes the wrong tool.** A single DP pairwise alignment is correct for two sequences. For one query against thousands to millions of targets (genome-scale homology search) or for many-vs-many all-by-all (clustering, ortholog detection), the right tool is profile- or k-mer-indexed search, not iterated DP:
- BLASTP / DIAMOND -- standard query-vs-database baseline
- MMseqs2 (Steinegger & Soding 2017 Nat Biotech) -- ~400x faster than PSI-BLAST at higher sensitivity; iterative profile mode `--num-iterations 3` matches PSI-BLAST and approaches `jackhmmer`
- MMseqs2-GPU (Kallenborn et al 2025 Nat Methods 22:2024; Mirdita co-author) -- GPU-accelerated; ~177x faster than `jackhmmer` for single queries on one NVIDIA L40S; use when GPU is available and the dataset is sensitivity-bound
- jackhmmer (HMMER) -- gold standard for distant homology when run to convergence; slow but the most sensitive non-structural method
- Foldseek -- escape to structural search when both query and database have predicted structures (see `alignment/structural-alignment`)

Other failure modes:
- **Non-homologous sequences**: All DP aligners return an alignment regardless of homology. E-value or bit score is the homology gate, not the existence of an alignment.
- **Repetitive sequences**: Tandem repeats produce ambiguous, artifactually high-scoring alignments; mask first.

## Related Skills

- alignment/multiple-alignment - Align three or more sequences with MAFFT, MUSCLE5, ClustalOmega
- alignment/alignment-io - Save alignments to files in various formats
- alignment/msa-parsing - Work with multiple sequence alignments
- alignment/msa-statistics - Calculate identity, similarity metrics
- alignment/structural-alignment - Twilight-zone alternative when sequence signal fails (Foldseek, TM-align, pLM aligners)
- alignment/alignment-trimming - Remove unreliable columns post-alignment
- sequence-manipulation/motif-search - Pattern matching in sequences

## References

- Karlin S, Altschul SF. 1990. Methods for assessing the statistical significance of molecular sequence features. PNAS 87:2264-2268.
- Schaffer AA et al. 2001. Improving the accuracy of PSI-BLAST protein database searches with composition-based statistics. NAR 29:2994-3005.
- Rost B. 1999. Twilight zone of protein sequence alignments. Prot Eng 12:85-94.
- Steinegger M, Soding J. 2017. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotech 35:1026-1028.
- Kallenborn F, Chacon A, Hundt C, Sirelkhatim H, Didi K, Cha S, Dallago C, Mirdita M, Schmidt B, Steinegger M. 2025. GPU-accelerated homology search with MMseqs2. Nat Methods 22(10):2024-2027.
- Sosic M, Sikic M. 2017. Edlib: a C/C++ library for fast, exact sequence alignment. Bioinf 33:1394-1395.
- Daily J. 2016. Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments. BMC Bioinf 17:81.
- Marco-Sola S et al. 2021. Fast gap-affine pairwise alignment using the wavefront algorithm. Bioinf 37:456-463.
- Marco-Sola S et al. 2023. Optimal gap-affine alignment in O(s) space. Bioinformatics 39(2):btad074.
- Yu YK, Altschul SF. 2005. The construction of amino acid substitution matrices for the comparison of proteins with non-standard compositions. Bioinf 21:902-911.
- Altschul SF, Erickson BW. 1985. Significance of nucleotide sequence alignments: a method for random sequence permutation that preserves dinucleotide and codon usage. MBE 2:526-538.
