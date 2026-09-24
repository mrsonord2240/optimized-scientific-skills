# Sequence Weighting and Neff

Read when the request involves sequence weights, redundancy correction or Neff. `find_conserved_positions` and `consensus_sequence` in SKILL.md accept `weights=`.

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
