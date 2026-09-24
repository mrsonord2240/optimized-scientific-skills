# Copy-Number Amplicon Bias Diagnostic

## Copy-Number Amplicon Bias Diagnostic

**Goal:** Detect the Aguirre 2016 / Munoz 2016 copy-number artifact where sgRNAs targeting amplified loci appear "essential" purely from DNA-damage burden.

**Approach:** Bin genes by copy number (if known from matched WGS/SNP-array) and check whether mean LFC correlates with CN. Alternatively, count off-target cut sites per sgRNA and check correlation with depletion -- amplified loci share many identical cut sites.

```bash
python scripts/cn_bias.py gene_lfc.tsv copy_number.tsv   # columns: gene,lfc and gene,copy_number
```

`scripts/cn_bias.py` (also importable: `cn_bias_diagnostic(gene_lfc_df, cn_df)`) merges the two tables, bins genes by copy number into quintiles, and returns the Spearman ρ of LFC vs copy number, the amplified (CN >4) vs diploid (CN 1.5-2.5) mean LFC and one-sided Mann-Whitney gap, `cn_bias_present` (either rule fires) and the per-bin table.

**Interpretation: two rules, not one.**

1. **Genome-wide.** Spearman ρ < -0.1 (p < 0.01) between copy number and LFC indicates a broad
   copy-number artifact.
2. **Focal.** Compare `amplified_mean_lfc` with `diploid_mean_lfc` directly. A single amplicon
   covers tens of genes out of ~18,000, so it barely moves ρ: on a realistic 40-gene amplicon the
   genome-wide ρ was only -0.066 while amplified genes averaged LFC -0.877 against -0.019 for
   diploid ones (p = 7e-19). Treat a gap below -0.5 with a significant one-sided test as bias even
   when ρ passes.

Either rule firing means correct before hit calling. When a specific amplicon is suspected, run the
diagnostic again on that region's genes plus a diploid background. Remediation: CRISPRcleanR, CERES
or Chronos (see [[copy-number-correction]], whose `detect_cn_bias()` applies the same two rules)
before hit calling.
