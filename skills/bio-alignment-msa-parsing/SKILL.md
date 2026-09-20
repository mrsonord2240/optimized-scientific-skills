---
name: bio-alignment-msa-parsing
description: Parse and analyze multiple sequence alignments using Biopython. Extract sequences, identify conserved regions, analyze gaps, work with annotations, and manipulate alignment data for downstream analysis. Use when parsing or manipulating multiple sequence alignments.
tool_type: python
primary_tool: Bio.AlignIO
license: MIT
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, numpy 1.26+, pyhmmer >= 0.11.3 (only for `compute_weights` and streaming). Checked 2026-09-19 on BioPython 1.88, numpy 2.0.2, pyhmmer 0.12.3.

Install: `pip install biopython numpy pyhmmer` (pyhmmer is optional).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# MSA Parsing and Analysis

Parse multiple sequence alignments to extract information, analyze content, and prepare for downstream analysis. Runnable versions of the helpers below are in `examples/` (`python examples/<name>.py [alignment_file]`; without an argument they run on the tiny alignment in `examples/data/`, and `msa_utils.py` holds the shared normalisation code).

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

HMMER, Stockholm and A2M files write gaps as `.` and soft-masked DNA is lowercase. Comparing against a literal `-` or counting letters case-sensitively then gives silently wrong answers (measured: gap counts all zero, `.` counted as residues, consensus `ACGTNNNN` on a soft-masked DNA alignment with unanimous columns). Every helper below normalises its input internally; call `normalize_alignment()` yourself before any raw `alignment[:, i]` comparison. For A2M/A3M pass `upper=False`, because case marks insert states there.

```python
def select_columns(alignment, keep, upper=None):
    '''New alignment with only the columns in `keep`; keeps record annotations,
    letter_annotations and column_annotations (so Stockholm GC/GR lines survive).
    upper=True/False also maps "." -> "-" (and upper-cases when True); None leaves sequences as is.'''
    keep = list(keep)
    records = []
    for record in alignment:
        seq = ''.join(str(record.seq)[i] for i in keep)
        if upper is not None:
            seq = (seq.upper() if upper else seq).replace('.', '-')
        new = SeqRecord(Seq(seq), id=record.id, name=record.name, description=record.description,
                        dbxrefs=list(record.dbxrefs), annotations=dict(record.annotations))
        for key, values in record.letter_annotations.items():
            picked = [values[i] for i in keep]
            new.letter_annotations[key] = ''.join(picked) if isinstance(values, str) else picked
        records.append(new)
    column_annotations = {}
    for key, values in getattr(alignment, 'column_annotations', {}).items():
        picked = [values[i] for i in keep]
        column_annotations[key] = ''.join(picked) if isinstance(values, str) else picked
    return MultipleSeqAlignment(records, annotations=dict(getattr(alignment, 'annotations', {})),
                                column_annotations=column_annotations)

def normalize_alignment(alignment, upper=True):
    return select_columns(alignment, range(alignment.get_alignment_length()), upper=upper)
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

Optional `weights` (one per sequence, e.g. from `henikoff_weights`) makes conservation phylogeny-aware; the denominator is then the sum of the weights.

```python
def find_conserved_positions(alignment, threshold=0.8, weights=None):
    alignment = normalize_alignment(alignment)
    weights = np.ones(len(alignment)) if weights is None else np.asarray(weights, dtype=float)
    if len(weights) != len(alignment):
        raise ValueError('weights must have one value per sequence')
    conserved = []
    for col_idx in range(alignment.get_alignment_length()):
        counts = Counter()
        for char, weight in zip(alignment[:, col_idx], weights):
            counts[char] += weight
        counts.pop('-', None)
        if not counts:
            continue
        most_common_char, most_common_count = counts.most_common(1)[0]
        conservation = most_common_count / weights.sum()
        if conservation >= threshold - 1e-12:
            conserved.append((col_idx, most_common_char, conservation))
    return conserved

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

## Alignment Trimming

Trimming controversy and tool selection (ClipKIT, trimAl, BMGE, Divvier, HMMcleaner, Noisy) is the subject of a dedicated skill. Use this short decision matrix for routing:

| Goal | First-line tool |
|------|-----------------|
| Phylogenetic-tree input | ClipKIT `kpic-smart-gap` (Steenwyk et al 2020 PLOS Bio) |
| HMM profile building | trimAl `-gappyout` (Capella-Gutierrez et al 2009 Bioinf) |
| Selection / dN/dS input | Avoid aggressive trimming; mask by MUSCLE5 ensemble column confidence (see Identifying Unreliable Alignment Regions) |
| Deep prokaryotic phylogenomics | BMGE (Criscuolo & Gribaldo 2010 BMC Evol Biol) |
| Preserve column-mapping for residue-level analysis | trimAl `-colnumbering` |

Prefer ClipKIT `kpic-smart-gap` over traditional gap-only removal for trees, and note that aggressive trimming (>20-30% of sites) can hurt tree quality. See alignment/alignment-trimming for full mode comparisons, decision trees, and runnable examples.

## Gap Handling for Phylogenetics

How gaps are treated in downstream phylogenetic analysis significantly affects tree topology:

| Treatment | Method | Tradeoff |
|-----------|--------|----------|
| Missing data (default) | Gaps = unknown character | Most common; can be statistically inconsistent under ML |
| Fifth state | Gap = 5th nucleotide | Biologically problematic (gaps of different lengths treated equally) |
| Simple indel coding | Each unique indel coded as binary character | Most biologically realistic; adds phylogenetic signal |

For slow- to mid-rate datasets where indels are phylogenetically informative, prefer SIC indel coding (use fifth-state only as a sensitivity check); for rapidly-evolving datasets (intra-species, ITS regions, retroelement-rich plant genomes), default to missing-data treatment because gap homology is unreliable. Run a sensitivity analysis comparing treatments before drawing topological conclusions.

## Identifying Unreliable Alignment Regions

Columns exhibiting **both** high gap fraction AND low conservation are the strongest indicators of alignment uncertainty. These often reflect guide tree artifacts rather than true evolutionary events; divergent sequences disproportionately introduce gaps. Before phylogenetic analysis:

1. Flag columns with gap fraction >50%, which may be alignment artifacts
2. Check if gappy regions coincide with insertions in a single divergent sequence (remove that sequence and re-align)
3. For critical analyses, get per-column confidence from a MUSCLE5 ensemble and mask low-confidence columns (checked on MUSCLE 5.3; needs the unaligned sequences, not the alignment). GUIDANCE2 is not offered here: its stand-alone package is no longer downloadable.

```bash
muscle -align seqs.fa -stratified -output ens.efa      # ensemble of replicate alignments
muscle -maxcc ens.efa -output maxcc.afa                # stderr ends "best <name>", e.g. acb.2
muscle -addconfseq ens.efa -output ens_cc.efa          # adds per-column confidence (CC) rows
python examples/muscle5_column_confidence.py ens_cc.efa acb.2 0.9 masked.fa
```

Run the script without the last two arguments first: it prints how many columns survive at CC 0.5/0.7/0.9/0.99, and there is no calibrated cut-off (GUIDANCE2's 0.93 does not transfer). On 8 UniProt globins 11 of 155 columns had CC < 0.9, including every column whose residue pairing differed in more than 10% of the 16 replicates (independent replicate-agreement check).

## Consensus Sequence

**"Get consensus sequence"** -> Derive a single representative sequence from an MSA based on majority-rule voting at each column.

**Goal:** Generate a consensus sequence from the alignment using a frequency threshold.

**Approach:** At each column, select the most common non-gap character if its share of all rows (gap rows included) reaches the threshold; otherwise mark as ambiguous. The placeholder is alphabet-aware: `N` for nucleotide, `X` for protein (`N` is asparagine and cannot be told apart from a real Asn column).

### Simple Majority Consensus
```python
def is_nucleotide(alignment, min_fraction=0.9):
    text = ''.join(str(r.seq) for r in alignment).upper().replace('-', '').replace('.', '')
    return bool(text) and sum(text.count(c) for c in 'ACGTUN') / len(text) >= min_fraction

def consensus_sequence(alignment, threshold=0.5, gap_char='-', ambiguous=None, weights=None):
    alignment = normalize_alignment(alignment)
    if ambiguous is None:
        ambiguous = 'N' if is_nucleotide(alignment) else 'X'
    weights = np.ones(len(alignment)) if weights is None else np.asarray(weights, dtype=float)
    if len(weights) != len(alignment):
        raise ValueError('weights must have one value per sequence')
    consensus = []
    for col_idx in range(alignment.get_alignment_length()):
        counts = Counter()
        for char, weight in zip(alignment[:, col_idx], weights):
            counts[char] += weight
        counts.pop('-', None)
        if not counts:
            consensus.append(gap_char)
            continue
        most_common_char, most_common_count = counts.most_common(1)[0]
        consensus.append(most_common_char if most_common_count / weights.sum() >= threshold else ambiguous)
    return ''.join(consensus)

consensus = consensus_sequence(alignment, threshold=0.5)
```

### Note on Bio.Align.AlignInfo
`AlignInfo.SummaryInfo` keeps only `get_column` in Biopython 1.88 (`dumb_consensus`, `gap_consensus`, `pos_specific_score_matrix`, `information_content` were removed and raise `AttributeError`). Use the custom `consensus_sequence()` above.

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

## Sequence Filtering

**Goal:** Subset an alignment to retain only sequences matching specific criteria (ID pattern, gap content, uniqueness).

**Approach:** Iterate over alignment records, apply filter conditions, and reconstruct a new MultipleSeqAlignment from matching records. A filter that would remove every sequence raises `ValueError` instead of returning an empty alignment.

```python
import re

def _require_kept(kept, alignment, what):
    if not kept:
        raise ValueError(f'{what} removes all {len(alignment)} sequences')
    return MultipleSeqAlignment(kept, annotations=alignment.annotations,
                                column_annotations=alignment.column_annotations)

def filter_by_id(alignment, pattern):
    regex = re.compile(pattern)
    return _require_kept([r for r in alignment if regex.search(r.id)], alignment, f'pattern {pattern!r}')

def filter_by_gap_content(alignment, max_gap_fraction=0.1):
    fractions = [str(r.seq).count('-') / len(r.seq) for r in normalize_alignment(alignment)]
    kept = [r for r, f in zip(alignment, fractions) if f <= max_gap_fraction]
    return _require_kept(kept, alignment, f'max_gap_fraction={max_gap_fraction} (lowest fraction {min(fractions):.2f})')

def remove_duplicates(alignment):
    seen = set()
    kept = [r for r in alignment if not (str(r.seq) in seen or seen.add(str(r.seq)))]
    return _require_kept(kept, alignment, 'remove_duplicates')
```

## Working with Annotations

Stockholm-derived alignments expose secondary-structure markup, GC/GR per-column annotations, and per-sequence metadata via `record.annotations`, `record.letter_annotations`, and `alignment.column_annotations`:

```python
alignment = AlignIO.read('pfam.sto', 'stockholm')
for record in alignment:
    if 'secondary_structure' in record.letter_annotations:
        print(record.id, record.letter_annotations['secondary_structure'])

ss_cons = alignment.column_annotations.get('secondary_structure')
```

GC SS_cons (consensus secondary structure), GC RF (reference coordinates), and GS metadata (organism, taxonomy) survive read/write through the `'stockholm'` format string but are silently discarded when writing to FASTA, PHYLIP, or NEXUS. Keep a Stockholm master copy if annotations matter for downstream analysis. `select_columns`, `remove_gappy_columns` and the filters above keep them; hand-built `SeqRecord(Seq(new_seq), id=...)` loops do not.

## Position Mapping

**Goal:** Convert between alignment column coordinates and ungapped sequence coordinates.

**Approach:** For one-off lookups, walk the sequence tracking gap characters. For repeated queries on the same sequence, vectorize with `numpy.cumsum` over a gap-mask -- O(L) preprocessing, O(1) lookups.

### Vectorized Coordinate Mapping (Recommended)

```python
import numpy as np

def coordinate_map(record):
    chars = np.frombuffer(str(record.seq).encode('ascii'), dtype=np.uint8)
    is_residue = ~np.isin(chars, [ord('-'), ord('.')])
    seq_to_aln = np.flatnonzero(is_residue)
    aln_to_seq = np.where(is_residue, np.cumsum(is_residue) - 1, -1)
    return seq_to_aln, aln_to_seq

seq_to_aln, aln_to_seq = coordinate_map(alignment[0])
column_of_residue_index_42 = seq_to_aln[42]   # 0-based index 42 is the 43rd residue
residue_index_at_column_100 = aln_to_seq[100]
```

`aln_to_seq[i] == -1` indicates a gap at alignment column `i`. This pattern handles 1 M-site genomic alignments in milliseconds compared to the loop-based version's seconds. It subsumes the single-lookup walk (`seq_pos += 1` for each non-gap char) in both directions.

### Mapping Alignment Columns to PDB Residues

A column-to-PDB mapping requires THREE coordinate systems: alignment column -> SEQRES residue (ungapped FASTA) -> ATOM residue (resolved structure). The SEQRES-to-ATOM map is non-trivial because PDB structures have unmodelled loops, N-terminal tags, engineered mutations, and seleno-substitutions. Conservation scores mapped via the bare alignment-to-SEQRES path will be off-by-many residues whenever the structure has missing density. Even the simplest case is offset: PDB 1MBN numbers residues from the first Val, so His93 is UniProt P02185 residue 94 (0-based index 93, checked). For SEQRES/ATOM extraction and the authoritative `_pdbx_poly_seq_scheme` mapping, see `structural-biology/structure-navigation`.

## Sequence Weighting and Neff

Compute sequence weights before column-wise statistics on phylogenetically structured datasets (lots of closely-related sequences plus a few outliers); without weighting, every per-column metric is biased toward the over-represented clades. `find_conserved_positions` and `consensus_sequence` accept `weights=`.

### Henikoff Sequence Weights

**Goal:** Give each sequence a weight that reflects its non-redundant contribution.

**Approach:** Each column `c` contributes `1 / (k_c * n_{s,c})` to sequence `s`, where `k_c` is the number of distinct residues at column `c` and `n_{s,c}` is the count of `s`'s residue at that column (Henikoff & Henikoff 1994 JMB). Columns containing any gap are skipped; weights sum to 1.

```python
def henikoff_weights(alignment):
    seq_array = np.array([list(str(r.seq).upper().replace('.', '-')) for r in alignment])
    weights = np.zeros(len(alignment))
    used_columns = 0
    for col_idx in range(seq_array.shape[1]):
        residues, inverse, counts = np.unique(seq_array[:, col_idx], return_inverse=True, return_counts=True)
        if '-' in residues:
            continue
        used_columns += 1
        weights += 1.0 / (len(residues) * counts[inverse])
    if not used_columns:
        raise ValueError('every column contains a gap, so no column can be weighted; '
                         'use pyhmmer compute_weights(method="pb") or trim gappy columns first')
    return weights / weights.sum()
```

Full implementation: `examples/henikoff_weights.py`.

**Edge case:** sequences whose residues fall ONLY in gap-containing columns receive weight zero (typically fragmentary or terminal-truncated sequences), and an alignment with no gap-free column raises `ValueError`. Easel's position-based weighting (`pyhmmer.easel` `compute_weights(method='pb')`, needs pyhmmer >= 0.11.3) does not skip columns: it ignores gaps column by column, considers only "consensus" columns (>= 50% residues by default), divides each weight by the sequence's residue count and rescales to sum to N. On the Pfam PF00042 seed it agrees with `henikoff_weights` in rank (Spearman 0.909) but not in value (max absolute difference 0.0091 after rescaling to sum 1). For gappy Pfam-style alignments use Easel:

```python
with pyhmmer.easel.MSAFile('alignment.sto', digital=True) as f:
    msa = f.read()
weights = np.array(msa.compute_weights(method='pb'), dtype=float)  # copy: pyhmmer reuses the vector
print(weights.sum())  # equals the number of sequences N, not an Neff
```

### Effective Sequence Number (Neff)

**Goal:** Estimate effective non-redundant sequence count after similarity-based clustering.

**Approach:** For each sequence count the sequences (itself included) whose identity over mutually non-gap positions is >= a threshold, give it weight `1 / count`, and sum. The protein default 0.62 is `hmmbuild --wblosum --wid`'s default (BLOSUM62's clustering threshold); 0.80 is a common DCA convention. Full implementation: `examples/neff.py`. This Skill's one rule of thumb for MI-APC and DCA contact prediction is **Neff/L > 1**; published cut-offs vary because each is calibrated against one estimator.

**Neff is estimator-dependent, and not every "weight" is an Neff.** Measured on the Pfam PF00042 seed (73 sequences, L = 141; HMMER 3.4, pyhmmer 0.12.3):

| Quantity | What it is | Value |
|----------|-----------|-------|
| `neff()` at identity 0.62 (`examples/neff.py`) | sum of 1/cluster-size: an Neff | 66.08 |
| pyhmmer `compute_weights('pb'/'gsc'/'blosum')` | per-sequence weights normalised to sum to N | 73.0 (sum), not an Neff |
| `hmmbuild` `eff_nseq` (default entropy weighting) | count tuned to a target relative entropy | 6.35 |

The first and last differ 10-fold on the same alignment. Always report which estimator was used, and apply a Neff/L threshold only to a number from the estimator it was calibrated on.

## Coevolution: Mutual Information with APC

**Goal:** Identify columns whose residue identities co-vary, indicating direct or indirect physical/functional coupling.

**Approach:** Compute pairwise mutual information across column pairs, then subtract the average-product correction (APC) from Dunn, Wahl & Gloor (2008 Bioinf) to remove per-column-entropy and phylogenetic background. This is the foundation of plmDCA and EVcouplings; full DCA needs Potts-model inference but APC-corrected MI runs in pure numpy and is informative on its own.

Skeleton (full implementation: `examples/mi_apc.py`):
```python
def mi_matrix_apc(alignment):
    # Compute pairwise MI across columns -> column_means -> APC = outer(means) / overall_mean
    # Return mi - apc
    ...
```

`examples/mi_apc.py` enforces the guard below itself: it computes Neff/L, warns and returns raw MI when the guard fails, and prints a column-shuffled null (each column permuted independently, best of 5 shuffles) so a ranking that does not beat it is visibly noise. Run it as `python examples/mi_apc.py alignment.sto`.

**Apply APC only when L > 100 and Neff/L > 1.** APC subtracts each column-pair's product of column-average MIs. For alignments with <100 columns, the column-averages are noisy estimates dominated by their constituent column-pairs, and APC removes signal proportional to noise; empirically on Pfam alignments <100 columns, APC-corrected MI underperforms raw MI for contact prediction (Cocco et al 2018 Rep Prog Phys review). Below the guard, raw MI plus a phylogenetic-distance threshold is more reliable, and shallow alignments give noise either way: on the PF00042 seed (L = 141, Neff/L = 0.47) the best MI-APC pair (0.603 bits) scored below the best column-shuffled pair (0.616), and none of the top 30 pairs were 1MBN contacts. For production-grade contact prediction, switch to plmDCA (Ekeberg et al 2013 Phys Rev E) or EVcouplings (Hopf et al 2017 Nat Biotechnol), which auto-skip APC when depth is insufficient. APC-corrected MI scales to a few hundred columns; deeper analyses need approximate likelihood methods.

## A2M / A3M Conventions

A2M (HMMER) and A3M (HHsuite, ColabFold) encode insert vs match columns via case (uppercase = match column residue, lowercase = insert). A3M does not pad inserts across sequences and must be reformatted to A2M before loading as a rectangular MSA. Match-only extraction: `examples/a2m_a3m_io.py`. See `alignment/alignment-io` A2M / A3M Conventions section for the full character table, BioPython load pattern, and the `reformat.pl` reference-sequence pitfall.

## Streaming Large Alignments

For Pfam-scale streaming (multi-gigabyte Stockholm or A3M databases that exceed RAM), use `pyhmmer.easel.MSAFile` with `compute_weights(method='pb')` for in-flight Henikoff weighting. See `alignment/alignment-io` Streaming Large Stockholm Databases section for the full code pattern.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `IndexError` | Column index out of range | Check `get_alignment_length()` |
| Unequal sequence lengths | Invalid MSA | Ensure all sequences same length |
| `ValueError: ... removes all N sequences` | Filter threshold too strict | Loosen the threshold; the message names the lowest gap fraction |
| `ValueError: every column contains a gap` | `henikoff_weights` on an all-gappy alignment | Trim gappy columns or use pyhmmer `compute_weights` |
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
