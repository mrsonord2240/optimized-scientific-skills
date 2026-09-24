---
name: bio-alignment-msa-parsing
category: Data Analysis
description: Parse and analyze multiple sequence alignments using Biopython. Extract sequences, identify conserved regions, analyze gaps, work with annotations, and manipulate alignment data for downstream analysis, including sequence weights, Neff, MI-APC coevolution and MUSCLE5 column confidence. Use when parsing or manipulating multiple sequence alignments.
tool_type: python
primary_tool: Bio.AlignIO
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, numpy 1.26+, pyhmmer >= 0.11.3 (only for `compute_weights` and streaming). Checked 2026-09-19 on BioPython 1.88, numpy 2.0.2, pyhmmer 0.12.3.

Install: `pip install biopython numpy pyhmmer` (pyhmmer is optional).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# MSA Parsing and Analysis

Parse multiple sequence alignments to extract information, analyze content, and prepare for downstream analysis. Runnable versions of the helpers below are in `examples/` (`python examples/<name>.py [alignment_file]`; without an argument they run on the tiny alignment in `examples/data/`, and `msa_utils.py` holds the shared normalisation code). Row filtering is `scripts/filter_sequences.py` (see `references/sequence-filtering.md`).

Conventions used throughout: column and residue positions are 0-based. **Conservation** is the fraction of sequences (gap rows included in the denominator) carrying the most common residue, so 80% in a 5-sequence alignment means less than 80% in a 500-sequence one; choose thresholds for the alignment's diversity.

## Required Import

**Goal:** Load modules for parsing, analyzing, and manipulating multiple sequence alignments.

**Approach:** Import AlignIO for reading, Counter for column analysis, and alignment classes for constructing modified alignments.

```python
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
import numpy as np
import pandas as pd
```

Optional for streaming and Easel-based weighting:
```python
import pyhmmer
```

## Loading Alignments

**Goal:** Read an MSA file and inspect its dimensions.

**Approach:** Use `AlignIO.read()` specifying the file and format.

```python
from Bio import AlignIO

alignment = AlignIO.read('alignment.fasta', 'fasta')
print(f'{len(alignment)} sequences, {alignment.get_alignment_length()} columns')  # alignment[i] is one record
```

## Gap and Case Normalisation

HMMER, Stockholm and A2M files write gaps as `.` and soft-masked DNA is lowercase. Comparing against a literal `-` or counting letters case-sensitively then gives silently wrong answers (measured: gap counts all zero, `.` counted as residues, consensus `ACGTNNNN` on a soft-masked DNA alignment with unanimous columns). Every helper below normalises its input internally; call `normalize_alignment()` yourself before any raw `alignment[:, i]` comparison. For A2M/A3M pass `upper=False`, because case marks insert states there (`references/a2m-a3m-streaming.md`).

`select_columns(alignment, keep, upper=None)` returns a new alignment holding only the column indices in `keep`, with record annotations, `letter_annotations` and `column_annotations` sliced to match (so Stockholm GC/GR lines survive); `upper=True/False` also maps `.` to `-` (and upper-cases when True), `None` leaves the sequences as they are. `normalize_alignment(alignment, upper=True)` is `select_columns` over every column. Both live in `examples/msa_utils.py`; the snippets below assume:

```python
import sys; sys.path.insert(0, 'examples')  # run from this Skill's directory
from msa_utils import normalize_alignment, select_columns
```

## Extracting Sequence Information

### Get All Sequence IDs
```python
seq_ids = [record.id for record in alignment]
```

### Get Sequences as Strings
```python
sequences = [str(record.seq) for record in alignment]
```

### Get Sequence by ID
```python
def get_sequence_by_id(alignment, seq_id):
    for record in alignment:
        if record.id == seq_id:
            return record
    return None

target = get_sequence_by_id(alignment, 'species_A')
```

### Access Descriptions and Annotations
```python
for record in alignment:
    print(f'ID: {record.id}')
    print(f'Description: {record.description}')
    print(f'Annotations: {record.annotations}')
```

## Column-wise Analysis

**Goal:** Analyze alignment content column by column to assess composition, conservation, and variability.

**Approach:** Use column indexing (`alignment[:, idx]`) and Counter to examine character frequencies at each position.

### Get Single Column
```python
column_5 = alignment[:, 5]  # Returns string of characters at position 5
print(column_5)  # e.g., 'AAAGA'
```

**API note:** `Bio.AlignIO` returns `MultipleSeqAlignment` objects whose `[:, idx]` returns a plain `str`; `[:, start:end]` returns another `MultipleSeqAlignment`. The newer `Bio.Align.Alignment` (from `Align.read` / `Align.parse`) uses numpy-backed slicing -- verify with `type(alignment[:, 0])` before assuming string methods work. For numpy-array access to the full alignment, use `np.array(alignment)`.

### Iterate and Count Columns
```python
for col_idx in range(alignment.get_alignment_length()):
    column = alignment[:, col_idx]
    counts = Counter(column)
```

### Find Conserved Positions

Optional `weights` (one per sequence, e.g. from `henikoff_weights`, see `references/weighting-neff.md`) makes conservation phylogeny-aware; the denominator is then the sum of the weights.

```python
from find_conserved import find_conserved_positions  # (column, residue, conservation) tuples; ValueError on wrong-length or all-zero weights

fully_conserved = find_conserved_positions(alignment, threshold=1.0)
mostly_conserved = find_conserved_positions(alignment, threshold=0.8)
```

## Gap Analysis

**Goal:** Quantify gap distribution across sequences and columns to identify problematic regions or sequences.

**Approach:** Count gap characters per sequence and per column, then identify positions exceeding a gap fraction threshold.

### Count Gaps Per Sequence
```python
gap_counts = [(record.id, str(record.seq).count('-')) for record in normalize_alignment(alignment)]
for seq_id, gaps in gap_counts:
    print(f'{seq_id}: {gaps} gaps')
```

### Count Gaps Per Column
```python
def gaps_per_column(alignment):
    alignment = normalize_alignment(alignment)
    return [alignment[:, i].count('-') for i in range(alignment.get_alignment_length())]

gap_profile = gaps_per_column(alignment)
```

### Find Gappy Columns
```python
def find_gappy_columns(alignment, threshold=0.5):
    profile = gaps_per_column(alignment)
    return [i for i, gaps in enumerate(profile) if gaps / len(alignment) >= threshold]

columns_to_remove = find_gappy_columns(alignment, threshold=0.5)
```

### Remove Gappy Columns
```python
def remove_gappy_columns(alignment, threshold=0.5):
    gappy = set(find_gappy_columns(alignment, threshold))
    keep = [i for i in range(alignment.get_alignment_length()) if i not in gappy]
    return select_columns(alignment, keep)  # annotations preserved

cleaned = remove_gappy_columns(alignment, threshold=0.5)
```

## Extracting Regions

### Slice by Column Range
```python
region = alignment[:, 100:200]  # Columns 100-199
```

### Slice by Sequence Range
```python
subset = alignment[0:10]  # First 10 sequences
```

### Extract Ungapped Regions from Reference
```python
def extract_ungapped_regions(alignment, ref_idx=0):
    ref_seq = str(normalize_alignment(alignment)[ref_idx].seq)
    return select_columns(alignment, [i for i, char in enumerate(ref_seq) if char != '-'])

ungapped = extract_ungapped_regions(alignment, ref_idx=0)
```

## Reference Files

Read the file when the request needs it; each holds runnable blocks moved verbatim from the sections named.

| File | Read when |
|------|-----------|
| `references/trimming-and-reliability.md` | The request is about trimming tools, how gaps enter a phylogenetic analysis, or masking unreliable columns (MUSCLE5 column confidence) (Alignment Trimming; Gap Handling for Phylogenetics; Identifying Unreliable Alignment Regions) |
| `references/consensus.md` | The request is for a consensus sequence (Consensus Sequence) |
| `references/sequence-filtering.md` | The request is to subset sequences by ID pattern, gap content or uniqueness (Sequence Filtering) |
| `references/annotations-and-position-mapping.md` | The request involves Stockholm annotations (SS_cons, RF, GS/GR) or converting between alignment columns and ungapped residue or PDB positions (Working with Annotations; Position Mapping) |
| `references/weighting-neff.md` | The request involves sequence weights, redundancy correction or Neff (Sequence Weighting and Neff) |
| `references/coevolution-mi-apc.md` | The request is to find coevolving column pairs (Coevolution: Mutual Information with APC) |
| `references/a2m-a3m-streaming.md` | The input is A2M/A3M or too large to load whole (A2M / A3M Conventions; Streaming Large Alignments) |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `IndexError` | Column index out of range | Check `get_alignment_length()` |
| Unequal sequence lengths | Invalid MSA | Ensure all sequences same length |
| `ValueError: ... removes all N sequences` | Filter threshold too strict | Loosen the threshold; the message names the lowest gap fraction |
| `ValueError: every column contains a gap` | `henikoff_weights` on an all-gappy alignment | Trim gappy columns or use pyhmmer `compute_weights` |
| `ValueError: weights must sum to a positive value` | all-zero `weights=` (e.g. from a degenerate weighting) | Recompute the weights; a 0/0 would otherwise return an all-placeholder consensus |
| Empty Counter | All gaps in column | Handled: consensus returns the gap character, conservation skips the column |

## Related Skills

- alignment/multiple-alignment - Run MSA tools (MAFFT, MUSCLE5, ClustalOmega) to generate alignments
- alignment/alignment-io - Read/write alignment files in various formats
- alignment/pairwise-alignment - Create pairwise alignments
- alignment/msa-statistics - Calculate conservation metrics
- alignment/alignment-trimming - ClipKIT, trimAl, BMGE, Divvier modes and decision trees
- alignment/structural-alignment - Twilight-zone alternative when sequence MSA is unreliable
- phylogenetics/modern-tree-inference - Build trees from processed alignments

## References

- Henikoff S, Henikoff JG. 1994. Position-based sequence weights. JMB 243:574-578.
- Dunn SD, Wahl LM, Gloor GB. 2008. Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction. Bioinf 24:333-340.
- Ekeberg M, Lovkvist C, Lan Y, Weigt M, Aurell E. 2013. Improved contact prediction in proteins: using pseudolikelihoods to infer Potts models. Phys Rev E 87:012707.
- Hopf TA, Ingraham JB, Poelwijk FJ, Scharfe CPI, Springer M, Sander C, Marks DS. 2017. Mutation effects predicted from sequence co-variation. Nat Biotechnol 35:128-135.
- Eddy SR. 2011. Accelerated profile HMM searches. PLOS CB 7:e1002195.
- Larralde M et al. 2023. PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinf 39:btad214.
- Cocco S et al. 2018. Inverse statistical physics of protein sequences: a key issues review. Rep Prog Phys 81:032601.
- Simmons MP, Ochoterena H. 2000. Gaps as characters in sequence-based phylogenetic analyses. Syst Biol 49:369-381.
- Mueller K. 2006. Incorporating information from length-mutational events into phylogenetic analysis. Mol Phylogenet Evol 38:667-676.
- Dwivedi B, Gadagkar SR. 2009. Phylogenetic inference under varying proportions of indel-induced alignment gaps. BMC Evol Biol 9:211.
- Velankar S et al. 2013. SIFTS: structure integration with function, taxonomy and sequences resource. NAR 41:D483-D489.
