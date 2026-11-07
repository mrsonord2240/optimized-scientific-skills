# Conditional or per-SNP methods

```bash
eCAVIAR -l ld_gwas.ld -l ld_eqtl.ld -z gwas.z -z eqtl.z -o coloc_out -c 2
```

Use ancestry-matched LD and aligned SNPs. Report maximum and sum CLPP as described in `references/ecaviar-clpp.md`. For conditional analysis with individual-level reference genotypes use `references/pwcoco.md`; report each conditional pair. SharePro requires a separate installed interface and validated input contract.

Done when the requested per-SNP or conditional results and assumptions are reported.
