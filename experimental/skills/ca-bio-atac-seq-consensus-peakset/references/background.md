# Background (not needed to run an analysis)

ATAC peak width varies across replicates because Tn5 cuts are stochastic at the edges: the same element can be 200 bp in one replicate and 800 bp in another. Counting reads in different-width intervals confounds width with signal, so differential analysis needs one fixed-width coordinate set that every sample is counted against.

Strategy comparison:

| Strategy | Width | Fails when |
|----------|-------|------------|
| Naive union (`bedtools merge`) | Variable, wide | Width inflation drives spurious differential |
| Intersection (all samples) | Variable | Loses condition-specific peaks |
| Majority rule (`minOverlap`) | Variable | Width still varies; unequal replicates drop peaks |
| Iterative overlap | 501 bp fixed | Loses sub-501 bp resolution; favours high-signal peaks |
| DiffBind summits | 501 bp at `summits=250` | Needs summits; broad peaks lose width |
| IDR-filtered union | Variable | Needs IDR per pair; mixed thresholds bias the set |
| Per-condition then union | Variable | Same width issues as a naive union |
| Midpoint extension | User-set | Midpoint is not the summit |

The strategy choice can shift peak counts and FDR by 30-50%, so document it. Verify against current ENCODE ATAC-seq standards before locking a pipeline. Iterative overlap yields a much smaller set than DiffBind summits because it drops overlapping summits; IDR yields fewer peaks than overlap counting because it is reproducibility-aware.

IDR thresholds differ: 0.05 for true replicates, 0.10 for pseudoreplicates.

Source reference material also showed a pybedtools re-centering function and a per-condition snippet; `scripts/iterative_overlap.sh`, `scripts/iterative_overlap.R` and `scripts/per_condition_consensus.sh` replace them.
