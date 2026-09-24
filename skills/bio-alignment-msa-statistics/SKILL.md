---
name: bio-alignment-msa-statistics
category: Data Analysis
description: Calculate alignment statistics including sequence identity, conservation scores, substitution matrices, and similarity metrics. Use when comparing alignment quality, measuring sequence divergence, and analyzing evolutionary patterns.
tool_type: python
primary_tool: Bio.AlignIO
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, numpy 1.26+ (checked on Biopython 1.88, numpy 2.0.2). Install: `pip install biopython numpy`.

Alignments must come from `AlignIO.read()` (`MultipleSeqAlignment`): the code below uses `alignment[:, i]` and `get_alignment_length()`, which `Bio.Align.read()` objects (`Alignment`) do not have.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# MSA Statistics

Calculate sequence identity, conservation scores, substitution counts, and other alignment metrics.

## Required Import and Normalisation

**Goal:** Load modules and put every alignment into one convention before computing anything.

**Approach:** Read with AlignIO, then normalise: `.`/`~` gaps become `-` and letters are upper-cased. MAFFT writes nucleotide alignments in lower case, and hmmalign/Stockholm/A2M use `.` gaps and lower-case inserts; every function below compares against `-` and counts letters case-sensitively, so un-normalised input gives **silently wrong numbers** (measured: Ti/Tv 0.00 instead of 1.21, DNA information content 29.9 bits instead of the 2.0-bit maximum, Pfam-seed mean PID1 32.0% instead of 19.6%). Run `check_alphabet()` afterwards to see which letters (X, B, Z, U, N) the alphabet-bound statistics will drop.

```python
import math
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, 'examples')   # this Skill's examples/ directory (run from the Skill directory, or use its absolute path)
from msa_utils import check_alphabet, load_alignment, normalize_alignment

alignment = load_alignment('alignment.fasta')   # AlignIO.read + normalize_alignment (U -> T for RNA); format guessed from the extension
check_alphabet(alignment, 'ACDEFGHIKLMNPQRSTVWY')  # 'ACGT' for DNA
```

`examples/msa_utils.py` holds `normalize_alignment` (`.`/`~` to `-`, upper-case, `u_to_t=True` for RNA; `upper=False` keeps case for A2M/A3M, where lower case marks insert columns; hmmalign A2M is ragged and cannot be loaded by AlignIO: see `alignment/alignment-io`), `check_alphabet` (prints the letters an alphabet-bound statistic will drop), `load_alignment`, `is_nucleotide`, `pick_background` and the backgrounds. `load_alignment` guesses the format from the extension (`.sto`/`.stk` Stockholm, `.aln`/`.clw` Clustal, `.phy` PHYLIP, `.nex` Nexus, otherwise FASTA) and takes `fmt=` to override, and the examples that read a file take that path as `argv[1]`. `is_nucleotide` calls an alignment DNA/RNA when at least 90% of its residues are A/C/G/T/U/N, or at least 50% are and the rest are IUPAC ambiguity codes (a DNA alignment with 12% R/Y/S/W/K/M is still DNA; protein sits near 30% A/C/G/T/N). Every example imports it and runs without arguments on the tiny alignments in `examples/data/`. `python examples/selftest.py` checks all statistics against values worked out by hand.

## Pairwise Identity

**"Calculate percent identity"** -> Compute the fraction of identical aligned residues between sequence pairs.

**Goal:** Measure sequence similarity as percent identity for individual pairs or across all sequences in an alignment.

**Approach:** Count matching non-gap positions divided by total aligned positions; optionally compute a full N-by-N identity matrix.

### Percent Identity Definitions

There are four common denominators, producing **up to 11.5% difference** on the same alignment. Combined with different alignment algorithms, variation reaches 22%. Always report which method was used.

| Method | Denominator |
|--------|-------------|
| PID1 | Aligned pairs + internal gap columns (the span from the pair's first to its last aligned column; terminal overhangs and unaligned flanks excluded, as in `pwalign::pid`) |
| PID2 | Aligned residue pairs only (no gaps) |
| PID3 | Shorter sequence length (ungapped) |
| PID4 | Mean sequence length (ungapped) |

PID2 always gives the highest value; PID4 correlates best with structural similarity (Raghava & Barton 2006 BMC Bioinf, r=0.86 with STAMP's Sc score) and is recommended for evolutionary analyses.

**Length-asymmetry pathology:** When sequences differ greatly in length, PID4 and PID2 diverge sharply. Example: 80 matches between a 500-residue protein and a 100-residue domain fragment yields PID2 ~84% (matches over aligned residue pairs) but PID4 ~27% (matches over mean ungapped length). Neither is wrong; they answer different questions:
- PID2 -> "how similar is the aligned region?" (motif/domain detection)
- PID4 -> "how similar are the full sequences?" (structural similarity benchmarks)

For ortholog identification at the protein level (full-length, similar size), PID4 is recommended. For domain detection or fragment-vs-genome alignment, PID2 with explicit length annotation is more interpretable. Always report alignment length alongside any percent identity to disambiguate.

### Calculate Identity Between Two Sequences
```python
from identity_matrix import pairwise_identity   # examples/identity_matrix.py; NaN when undefined (all-gap sequence, no aligned pair for pid1/pid2)

seq1, seq2 = str(alignment[0].seq), str(alignment[1].seq)
for method in ['pid1', 'pid2', 'pid3', 'pid4']:
    print(f'{method}: {pairwise_identity(seq1, seq2, method) * 100:.1f}%')
```

An X/X or N/N column counts as identical. Cross-check PID1 with `pwalign::pid(aln, type='PID1')` (R, Bioconductor) when a number matters.

### Identity Matrix for All Sequences

The double-loop is O(N^2 * L) and fine for hundreds of sequences; for thousands, vectorize via numpy broadcasting:

```python
def identity_matrix_vectorized(alignment, method='pid4'):
    # Build N x L character array; for each row, broadcast equality and validity masks against all rows; NaN where undefined
    ...
```

Full implementation (all four PID methods, default `pid4`, NaN-aware `average_identity`): `examples/identity_matrix.py [alignment] [pid1..pid4]`. For very large alignments (>10k sequences), switch to k-mer-based distance estimation (e.g. mash) -- exact pairwise identity becomes prohibitive.

## Conservation Scoring Methods

Pick a conservation score by what the downstream task needs:

| Pick this | When |
|-----------|------|
| Majority fraction or Shannon entropy | Quick screening; DNA/RNA logos; coarse column QC (entropy: `references/information-content-pssm.md`) |
| Capra-Singh JSD (modern default) | Catalytic-residue / functional-site prediction on protein MSAs (`references/capra-singh-jsd.md`) |
| ConSurf rate4site | PDB surface mapping when a phylogenetic tree is available |

## Conservation Score

**Goal:** Quantify per-column and overall alignment conservation to identify conserved and variable regions.

**Approach:** Calculate the fraction of the most common residue at each column, ignoring gaps, and smooth with a centred sliding window. Occupancy matters: ignoring gaps makes a column with 2 residues out of 73 sequences "100% conserved", so columns with fewer than `min_occupancy` (default 0.5) residues, and all-gap columns, score NaN and are skipped by the averages (on the Pfam globin seed this moves the mean from 37.9% to 34.3%).

### Per-Column Conservation
```python
def column_conservation(alignment, col_idx, ignore_gaps=True, min_occupancy=0.5):
    full = alignment[:, col_idx]
    column = full.replace('-', '') if ignore_gaps else full
    if not column or (ignore_gaps and len(column) / len(full) < min_occupancy):
        return float('nan')
    return Counter(column).most_common(1)[0][1] / len(column)

for i in range(min(20, alignment.get_alignment_length())):
    cons = column_conservation(alignment, i)
    print(f'Column {i}: {cons*100:.0f}% conserved')
```

**Ranking columns:** NaN columns make `sorted(..., key=lambda i: -scores[i])` and `max()` silently misorder (NaN compares False with everything; measured on the Pfam globin seed, the "top 10" then holds values as low as 0.37 instead of 0.63-1.0). Filter NaN first:

```python
scores = [column_conservation(alignment, i) for i in range(alignment.get_alignment_length())]
top10 = sorted((i for i, s in enumerate(scores) if not math.isnan(s)), key=lambda i: -scores[i])[:10]
```

### Average Conservation Across Alignment
```python
from conservation_profile import average_conservation   # examples/conservation_profile.py; (nan, 0) when no column reaches min_occupancy

avg_cons, n_used = average_conservation(alignment)
if n_used:
    print(f'Average conservation: {avg_cons*100:.1f}% over {n_used} columns')
else:
    print('No column has >= min_occupancy residues: lower min_occupancy or trim the alignment first')
```

### Conservation Profile
```python
def conservation_profile(alignment, window=10, min_occupancy=0.5):
    # centred moving average over columns i-window//2 .. i+window//2 (inclusive), NaN columns skipped
    length = alignment.get_alignment_length()
    scores = np.array([column_conservation(alignment, i, min_occupancy=min_occupancy) for i in range(length)])
    profile = []
    for i in range(length):
        chunk = scores[max(0, i - window // 2):i + window // 2 + 1]
        chunk = chunk[~np.isnan(chunk)]
        profile.append(chunk.mean() if chunk.size else float('nan'))
    return profile

profile = conservation_profile(alignment, window=10)
```

### Capra-Singh Jensen-Shannon Divergence

JSD conservation with window smoothing for catalytic-residue prediction (Capra & Singh 2007): `references/capra-singh-jsd.md`. Implementation: `examples/capra_singh_jsd.py`.

## Substitution Counts

Raw pairwise substitution counts and Ti/Tv from an alignment, `PairwiseAligner.substitutions`, and why the BLOSUM62 lambda is not one number: `references/substitution-counts.md`. Implementation: `examples/substitution_counts.py`.

## Information Content, PSSM, Neff, MI-APC

Per-column Shannon entropy, information content as KL divergence against the Robinson (protein) or uniform (DNA) background, pseudocount PSSM, Neff and MI-APC: `references/information-content-pssm.md`. Implementations: `examples/entropy_analysis.py`, `examples/pssm.py`. `ROBINSON_BACKGROUND` and `DNA_UNIFORM` are defined in `examples/msa_utils.py`.

## Gap Statistics

Gap fraction per column and a summary (total gaps, gappiest sequence and column, gap-free columns): `references/gap-statistics.md`. Implementation: `examples/gap_statistics.py`.

## Alignment Quality Metrics

Flat-gap `alignment_score` and BLOSUM62 `sum_of_pairs` (gap/gap pairs score 0; SP bias on unbalanced datasets): `references/alignment-quality-scores.md`. Implementation: `examples/alignment_scores.py`.

## Distance Correction Models

For publication-grade distances select the model with ModelTest-NG, then correct with IQ-TREE2 or `distmat`; `DistanceCalculator` is exploratory only: `references/distance-correction.md`.

## Alignment Quality Assessment

### When to Worry About Alignment Quality

| Warning Sign | Implication | Action |
|-------------|-------------|--------|
| Average pairwise identity <25% (protein) | Twilight zone; alignment may be unreliable | Quantify per-column confidence (see below); consider structural alignment |
| >30% of columns have >50% gaps | Possible non-homologous sequences, misalignment or guide-tree artifacts | Remove outlier sequences and re-align |
| Identity varies dramatically across regions | Domain architecture mismatch | Align domains separately |
| Conservation pattern absent in expected functional regions | Alignment error or non-homology | Verify with BLAST that sequences are truly homologous |

### Quantifying Alignment Uncertainty

Alignment uncertainty propagates directly into evolutionary inference; for selection analysis (dN/dS), misaligned codons create artificial nonsynonymous differences and false positive signals. Different aligners can support different tree topologies: report the alignment method and run a sensitivity analysis. For critical analyses (phylogenetics, selection) quantify per-column confidence first (GUIDANCE2, MUSCLE5 ensemble, T-Coffee TCS): see the Confidence Assessment section in `alignment/multiple-alignment`.

## Note on Bio.Align.AlignInfo

The `AlignInfo.SummaryInfo` class is **deprecated** in recent Biopython versions. Use the custom functions in this skill instead:
- For PSSM: use `pssm_with_pseudocounts()` (`references/information-content-pssm.md`)
- For information content: use the `information_content()` function (`references/information-content-pssm.md`)
- For consensus: see msa-parsing skill

## Quick Reference: Metrics

| Metric | Description | Range | Higher means |
|--------|-------------|-------|--------------|
| Identity | Fraction of identical residues (choose PID1-4) | 0-1 | more similar |
| Conservation | Most common residue frequency | 0-1 | less variable |
| Shannon Entropy | Variability measure | 0 to log2(alphabet) | more variable |
| Information Content | KL divergence from background (equals max entropy - observed only for a uniform background) | 0 to -log2(rarest background letter): 2 bits DNA, ~6.2 bits protein | more constrained |
| Gap Fraction | Proportion of gaps | 0-1 | less reliable column |

## Common Errors

| Error / symptom | Cause | Solution |
|-------|-------|----------|
| Ti/Tv 0.00, IC far above its maximum, identity too low, no error raised | Lower-case letters or `.` gaps not normalised | `normalize_alignment()` first; `check_alphabet()` |
| `IndexError` in `substitution_matrix[c1, c2]` | Letter outside BLOSUM62 (`U`, `J`, lower case, `.`) | Normalise; `sum_of_pairs` skips and counts them |
| `ValueError: Bad letter` from `DistanceCalculator` | Lower case, `.` or ambiguity letters | Normalise; drop or mask ambiguity codes first |
| `AttributeError: 'Alignment' object has no attribute 'get_alignment_length'` | Alignment read with `Bio.Align.read()` | Read with `AlignIO.read()` |
| NaN identity or conservation | Undefined by definition (all-gap sequence, column below `min_occupancy`) | Expected; averages skip NaN |

## Reference Files

Read the file that matches the request; each assumes the normalised `alignment` from the section above.

| File | Read when |
|------|-----------|
| `references/capra-singh-jsd.md` | JSD conservation / catalytic-residue prediction |
| `references/substitution-counts.md` | Substitution pair counts, Ti/Tv, `.substitutions`, BLOSUM62 lambda |
| `references/information-content-pssm.md` | Shannon entropy, information content (KL), PSSM, Neff, MI-APC |
| `references/gap-statistics.md` | Gap fraction per column, gap summary |
| `references/alignment-quality-scores.md` | `alignment_score`, BLOSUM62 sum of pairs, SP bias |
| `references/distance-correction.md` | ModelTest-NG / IQ-TREE corrected distances |
| `references/bibliography.md` | Citing a definition, background or threshold |

## Related Skills

- alignment/multiple-alignment - Run MSA tools and quantify alignment confidence with ensembles
- alignment/msa-parsing - Parse, filter, trim, and assess alignment quality (Henikoff weights, Neff, MI-APC live there)
- alignment/alignment-io - Read/write alignment files
- alignment/pairwise-alignment - Create and score pairwise alignments
- alignment/alignment-trimming - Column trimming before downstream statistics
- alignment/structural-alignment - Twilight-zone alternative when sequence MSA is unreliable
- phylogenetics/distance-calculations - Distance models and tree building from corrected distances
- sequence-manipulation/sequence-properties - Sequence-level statistics
