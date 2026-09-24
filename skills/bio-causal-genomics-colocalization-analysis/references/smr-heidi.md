# SMR + HEIDI

## SMR vs coloc Reconciliation

SMR (Zhu 2016) and coloc test related but non-identical questions:

- **SMR** tests pleiotropy vs linkage: does the top eQTL SNP show a GWAS effect explainable by its eQTL effect (pleiotropic / causal) or does the GWAS effect come from a different SNP in LD (linkage)?
- **coloc** tests shared vs distinct causal variants over an entire window of SNPs.
- **HEIDI** is SMR's heterogeneity test; null hypothesis is single shared causal SNP. Zhu 2016 Nat Genet 48:481 specifies **HEIDI p > 0.05** (NOT 0.01) as non-rejection of single shared causal. HEIDI p > 0.05 does NOT prove shared causality -- only that data cannot reject it; pair with SMR p Bonferroni-corrected across probes. When LD between causal SNPs > 0.7, HEIDI loses power same as coloc.

When LD between two true causal SNPs is high (r2 > 0.7), both SMR/HEIDI and coloc.abf lose discriminatory power: SMR cannot pick which of the LD-tied SNPs is causal, and coloc.abf cannot reject H4 even if biology is two-distinct-causal. coloc.susie + ancestry-matched LD is the modern resolution.

**Operational rule:** SMR + HEIDI is appropriate when the question is "does this eQTL gene mediate the GWAS effect at all?" coloc is appropriate when the question is "do the two traits share a causal variant in this window?" Run both; agreement (significant SMR + non-rejected HEIDI + PP.H4 >= 0.75) is high-confidence; disagreement requires inspection (often the multi-causal / LD scenario above).

## SMR + HEIDI Pipeline

```bash
# SMR is a command-line tool. Pre-format GWAS into .ma (SNP A1 A2 freq beta se p N).
# eQTL data as BESD (binary eQTL summary data); pre-built BESD available from eQTLGen / GTEx.

smr --bfile 1KG_EUR_chr6 \
    --gwas-summary gwas.ma \
    --beqtl-summary eqtl_chr6.besd \
    --out smr_result \
    --thread-num 4 \
    --peqtl-smr 5e-8 \
    --heidi-mtd 1
# Output smr_result.smr: probe (gene) | top SNP | p_SMR | p_HEIDI | nsnp_HEIDI
```

Interpretation: significant `p_SMR` (Bonferroni-corrected across probes tested, typically < 5e-8 / N_probes) AND non-rejection by HEIDI (`p_HEIDI > 0.05`, per Zhu 2016) indicates pleiotropy / shared causal; `p_HEIDI <= 0.05` rejects shared-causal -> linkage. HEIDI p > 0.05 does NOT prove shared causality, only that data cannot reject it. Require `nsnp_HEIDI >= 10` for HEIDI reliability.
