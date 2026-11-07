# SMR and HEIDI

```bash
smr --bfile reference --gwas-summary gwas.ma --beqtl-summary eqtl.besd --out smr_result --thread-num 4 --peqtl-smr 5e-8 --heidi-mtd 1
```

See `references/smr-heidi.md` for inputs and interpretation. Correct SMR across probes, report HEIDI SNP count and require at least 10 for its reliability; HEIDI p > 0.05 is non-rejection, not proof of shared causality.

Done when SMR, HEIDI and limitations are reported.
