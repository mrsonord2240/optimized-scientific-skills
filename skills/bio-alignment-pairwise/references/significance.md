# Statistical Significance (reference)

## Statistical Significance: Karlin-Altschul

Use bit score (database-size-independent) and E-value (expected chance hits) to interpret raw alignment scores; never compare raw scores across scoring schemes. For non-default gap penalties or non-protein/DNA alphabets, generate empirical p-values via sequence shuffling instead of trusting the formula (`examples/empirical_pvalue.py`; for DNA pass `preserve='di'`, the dinucleotide shuffle of Altschul & Erickson 1985, pure Python, no `ushuffle` build).

| Bit score | E-value (typical 1e6 db) | Interpretation |
|-----------|--------------------------|----------------|
| > 200 | < 1e-50 | Essentially certain homology |
| 50-200 | 1e-5 to 1e-50 | Likely homology |
| 30-50 | 0.01 to 1e-5 | Possible homology; verify with profile methods |
| < 30 | > 0.01 | Suspect; not significant |

Karlin-Altschul lambda/K are calibrated empirically per (matrix, gap-open, gap-extend) tuple; switching gap penalties without recalibrating produces wrong E-values.

**Mask compositional bias before scoring.** Low-complexity regions (poly-Q, proline-rich, leucine-rich TM helices) inflate raw scores. Pre-filter with SEG for protein or DUST for DNA, then run BLAST with `-comp_based_stats 2` (modern default; Yu & Altschul 2005 Bioinf, conditional) or `-comp_based_stats 3` for repeat-rich queries (Yu & Altschul 2005 unconditional; Schaffer et al 2001 NAR introduced the original level-1 statistics). SEG pre-filtering is complementary to `-comp_based_stats`, not redundant.

## References

- Karlin S, Altschul SF. 1990. Methods for assessing the statistical significance of molecular sequence features. PNAS 87:2264-2268.
- Schaffer AA et al. 2001. Improving the accuracy of PSI-BLAST protein database searches with composition-based statistics. NAR 29:2994-3005.
- Yu YK, Altschul SF. 2005. The construction of amino acid substitution matrices for the comparison of proteins with non-standard compositions. Bioinf 21:902-911.
- Altschul SF, Erickson BW. 1985. Significance of nucleotide sequence alignments: a method for random sequence permutation that preserves dinucleotide and codon usage. MBE 2:526-538.
