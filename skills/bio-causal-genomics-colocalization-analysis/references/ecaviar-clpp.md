# eCAVIAR CLPP

## eCAVIAR CLPP Threshold Framework

CLPP (Colocalization Posterior Probability) is the per-SNP product of the two per-trait fine-mapping posteriors. Threshold conventions:

- Hormozdiari 2016 AJHG 99:1245 used CLPP >= 0.01 (validated against null simulations).
- 2024 GTEx / Open Targets pipelines use CLPP >= 0.05.
- High-confidence claims require CLPP >= 0.1.
- Report both sum-CLPP across the credible set AND max-CLPP at any single SNP -- the two answer different questions (locus-level vs lead-SNP-level confidence).

```bash
eCAVIAR -l ld_gwas.ld -l ld_eqtl.ld \
        -z gwas.z -z eqtl.z \
        -o coloc_out -c 2     # -c = max independent causal variants per trait
# Output: per-SNP CLPP in coloc_out_col file; report sum and max
```
