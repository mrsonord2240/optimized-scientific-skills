# Information content, PSSM, Neff and MI-APC

Read when computing Shannon entropy, information content (KL against a background), a pseudocount PSSM, Neff or coupling (MI-APC). Assumes the normalised `alignment` and the imports from "Required Import and Normalisation" in `SKILL.md`.

## Information Content

**Goal:** Measure column variability using Shannon entropy and derive information content for identifying functionally important positions.

**Approach:** Compute Shannon entropy from character frequencies per column; information content is the divergence of the column from a background (below). Letters outside the background (X, B, Z, U, N) are dropped and the column renormalised: giving them a tiny fallback probability inflates IC by ~26 bits per unknown letter (measured 29.9 bits on a DNA alignment whose maximum is 2.0).

### Shannon Entropy Per Column
```python
from entropy_analysis import shannon_entropy   # examples/entropy_analysis.py

for i in range(min(20, alignment.get_alignment_length())):
    column = alignment[:, i]
    ent = shannon_entropy(column)
    print(f'Column {i}: entropy = {ent:.2f} bits')
```

### Information Content (Kullback-Leibler Divergence)

The classic uniform-background formulation `IC = log2(alphabet_size) - H` (Schneider & Stephens 1990 NAR) is only valid when the genomic background is uniform. This is approximately true for random DNA but emphatically wrong for protein, where amino acid frequencies range from 1.3% (Trp) to 9.0% (Leu). For amino acids, use Kullback-Leibler divergence `IC = sum_i p_i * log2(p_i / b_i)` against the Robinson & Robinson 1991 PNAS empirical background (NCBI-tabulated values in `examples/msa_utils.py`; sum = 1.0). Full implementation: `examples/entropy_analysis.py`, which chooses the background from the alphabet.

```python
from entropy_analysis import information_content   # examples/entropy_analysis.py: drops gaps and unknown letters, renormalises the rest
from msa_utils import pick_background   # examples/msa_utils.py: DNA_UNIFORM for nucleotide, ROBINSON_BACKGROUND for protein

background, label = pick_background(alignment)
print([round(information_content(alignment[:, i], background), 2) for i in range(10)])
```

`python examples/entropy_analysis.py [alignment]` prints entropy and IC per column and chooses the background from the alphabet.

For sequence-logo letter heights, use the Schneider-Stephens form (`letter_height = p_i * (log2(alphabet) - H_observed)`, uniform background) when the comparison is "informative vs random"; use the KL form for protein logos or when the comparison is "informative vs the proteome". When the background is unknown, default to the empirical alignment composition rather than uniform.

## Position-Specific Score Matrix (PSSM)

**Goal:** Build a position-specific scoring matrix from the alignment for motif analysis or sequence scoring.

**Approach:** Raw counts give frequencies; without pseudocounts, log-odds against background diverge to negative infinity at any column missing a residue. Henikoff JG & Henikoff S 1996 (Bioinf 12:135-143) introduced data-dependent pseudocount weighting; a simple total pseudocount of 1 spread over residues by background frequency (not a per-residue Laplace add-one) is the minimal correct approach for production use. Letters outside the background are dropped from the counts and the column total. Full implementation: `examples/pssm.py`, which picks a protein or nucleotide background from the alignment (the A/C/G/T letters exist in both alphabets, so a protein background would otherwise run silently on DNA).

```python
def pssm_with_pseudocounts(alignment, background, pseudocount=1.0):
    # log2((counts[r] + pc * background[r]) / (n + pc) / background[r]) per column, n = residues in the background
    ...
```

`pseudocount=1.0` is the total pseudocount per column; HMMER uses Dirichlet mixtures for sophisticated smoothing. For motif scanning, score a candidate site by summing per-position log-odds; sites above a calibrated threshold are predicted hits. Use `ROBINSON_BACKGROUND` (`examples/msa_utils.py`) for protein.

## Effective Sequence Number (Neff)

**Goal:** Estimate non-redundant sequence count for MSA-depth metrics.

**Approach:** Cluster at an identity threshold (0.62 protein, 0.80 nucleotide) and weight by inverse cluster size; reference implementation lives in `msa-parsing` (`examples/neff.py`). `Neff/L > 0.5` is the rule-of-thumb for direct-coupling-analysis contact prediction; AlphaFold's MSA-depth scoring uses a closely related metric.

## Mutual Information with APC

**Goal:** Detect coevolving column pairs as a coupling/contact signal.

**Approach:** Pairwise MI minus average-product correction (Dunn, Wahl, Gloor 2008 Bioinf). Reference implementation lives in `msa-parsing` (`examples/mi_apc.py`). For production-grade contact prediction beyond a few hundred columns, switch to plmDCA (Ekeberg et al 2013) or EVcouplings (Hopf et al 2017).
