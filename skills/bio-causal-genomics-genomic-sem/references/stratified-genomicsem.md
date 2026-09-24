# Stratified GenomicSEM (reference for genomic-sem SKILL.md)

Partitioned heritability of the common factor with `enrich()`. Read it when `SKILL.md` "Reference Files" points here. The section names below refer to `SKILL.md`.

## Stratified GenomicSEM (Partitioned Heritability of Factor)

For partitioning the heritability of the latent factor across functional annotations, use `s_ldsc()` (stratified LDSC inside GenomicSEM) and pass the multi-annotation output to a stratified model fit.

```r
# Stratified LDSC across baseline + custom annotations
s_results <- s_ldsc(
    traits = traits,
    sample.prev = c(0.5, 0.5, NA),
    population.prev = c(0.05, 0.05, NA),
    ld = 'baselineLD_v2.2.',
    wld = 'weights.hm3_noMHC.',
    frq = '1000G.EUR.QC.',
    trait.names = trait_names
)

# enrich() inventory:
#   params: lavaan syntax of the parameter under enrichment (loading, residual var, or F~~F latent var)
#   fix='regressions': hold regression paths fixed at the genome-wide estimate during stratified fit
#   std.lv=FALSE: do not standardize the latent variance
#   rm_flank=TRUE: drop flanking-window contributions (default)
#   tau=FALSE: use the baseline-annotation S/V matrices (TRUE switches to the V_Tau/S_Tau tau parametrization)
#   base=TRUE: include baseline annotation contributions in the partition
#   toler=NULL: matrix-inversion tolerance (let GenomicSEM choose; supply a small value when S is near-singular)
strat_factor <- enrich(s_covstruc = s_results,
                       model = '',
                       params = 'F =~ trait1',
                       fix = 'regressions',
                       std.lv = FALSE,
                       rm_flank = TRUE,
                       tau = FALSE,
                       base = TRUE,
                       toler = NULL)
```

The output gives per-annotation enrichment of the factor h2 -- the analog of cell-type S-LDSC for the latent factor (Grotzinger AD et al 2022 Nat Genet 54:548).
