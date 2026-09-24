# Coevolution: Mutual Information with APC

Read when the request is to find coevolving column pairs. Check Neff first (weighting-neff.md).

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

**Apply APC only when L > 100 and Neff/L > 1.** APC subtracts each column-pair's product of column-average MIs. For alignments with <100 columns, the column-averages are noisy estimates dominated by their constituent column-pairs, and APC removes signal proportional to noise; treat APC as unreliable there (this is reasoning from the estimator, not a measured, cited result). Below the guard, raw MI plus a phylogenetic-distance threshold is more reliable, and shallow alignments give noise either way: on the PF00042 seed (L = 141, Neff/L = 0.47) the best MI-APC pair (0.603 bits) scored below the best column-shuffled pair (0.616), and none of the top 30 pairs were 1MBN contacts (top-30 contact precision: raw MI 0.13, MI-APC 0.00, random baseline 0.09, so both are noise). For production-grade contact prediction, switch to plmDCA (Ekeberg et al 2013 Phys Rev E) or EVcouplings (Hopf et al 2017 Nat Biotechnol), which auto-skip APC when depth is insufficient. Background on APC and DCA: Cocco et al 2018 (review). APC-corrected MI scales to a few hundred columns; deeper analyses need approximate likelihood methods.
