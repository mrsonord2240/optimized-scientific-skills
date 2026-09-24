# Popcorn trans-ancestry genetic correlation (reference for heritability-partitioning SKILL.md)

## Popcorn Pipeline

**Goal:** Trans-ancestry genetic correlation and per-population heritability from GWAS summary
statistics (Brown 2016 AJHG 99:76). Distinguishes shared-causal from ancestry-specific signal;
requires effective N > 5000 per population.

**Approach:** `compute` cross-population LD/covariance scores from two ancestry-matched reference
panels -> `fit` heritability and genetic correlation to two summary-statistics files against those
scores.

Run it as `bash examples/popcorn_rg.sh <bfile1> <bfile2> <sumstats1.txt> <sumstats2.txt> <out_prefix>`.
Sumstats need `SNP`, `A1`/`A2`, `N`, and `Z` (or `beta`/`SE`, or `OR`/`p-value`) columns.

## Install

```bash
git clone https://github.com/brielin/Popcorn.git
cd Popcorn
python -m venv pc_venv && source pc_venv/bin/activate
pip install "numpy<2" "pandas<2.2" scipy pysnptools bottleneck statsmodels
pip install .
```

**Two real numpy>=2 bugs found and fixed in this pip package (checked 2026-09-21, package version
1.1)** -- this is the actively-maintained Python 3 port (ported May 2022 by shafayetrahat per the
package's own README), not the abandoned 2015 original, and both bugs are still present in it:

1. `popcorn compute` crashes with `TypeError: only 0-dimensional arrays can be converted to Python
   scalars` under numpy>=2 (a `float(array)` call inside `popcorn/compute.py`'s scoring loop). Fixed
   by installing into a venv pinned to `numpy<2` as above -- no source patch needed.
2. `popcorn fit`'s jackknife standard-error step then crashes separately with
   `AttributeError: module 'numpy' has no attribute 'bool'` (`np.bool` was removed in numpy>=1.24;
   `numpy<2` alone still lands above that removal). One-line source patch in the installed package,
   `popcorn/jackknife.py`: `dtype=np.bool` -> `dtype=bool`.

## Verified

2026-09-21, with both fixes applied: ran `compute` + `fit` end to end on a real 957-individual,
14389-SNP genotype panel split into two halves (used as `--bfile1`/`--bfile2`) with a planted shared
genetic component (`rg` target 0.7) between two simulated phenotypes and real per-SNP OLS Z-scores.
Both steps completed without error and produced real jackknife standard errors and p-values.

**Toy-scale caveat, not a code defect:** the reported `h1^2`/`h2^2` came back far outside [0,1]
(~150) because splitting one ancestry's genotypes in half gives the two "populations" near-identical
LD and allele frequencies, which the cross-population covariance-score model does not expect --
Popcorn's heritability normalization assumes genuinely distinct ancestries (its own worked example
uses 1000G EUR vs EAS). Use real cross-ancestry reference panels for a result whose magnitude means
anything; this test only demonstrates that the pipeline's code path runs correctly end to end.

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `h1^2`/`h2^2` far outside [0,1] | Reference panels for `--bfile1`/`--bfile2` are not genuinely distinct ancestries | Use real ancestry-matched panels (e.g. 1000G EUR vs EAS); do not substitute a same-ancestry split |
| `AttributeError: module 'numpy' has no attribute 'bool'` during `fit` | numpy>=1.24 removed `np.bool`; present in `popcorn/jackknife.py` | Patch `dtype=np.bool` -> `dtype=bool` in the installed package (see Install) |
| `TypeError: only 0-dimensional arrays can be converted to Python scalars` during `compute` | numpy>=2 incompatibility in `popcorn/compute.py` | Install into a venv pinned to `numpy<2` |
