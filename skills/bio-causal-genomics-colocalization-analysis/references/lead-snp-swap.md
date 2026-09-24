# Window-centring diagnostics

### Lead-SNP swap and window bias

**Trigger:** The two traits have different lead SNPs at the same locus; analyst centres each window on the trait-specific lead.

**Mechanism:** coloc PP is sensitive to the SNPs in the window; centring on different leads gives different per-SNP overlap and biases toward H3.

**Symptom:** Re-centring the window on the GWAS lead vs the eQTL lead produces qualitatively different PP.H4.

**Fix:** Use a SINGLE window (typically +/- 500 kb or 1 Mb) centred on the joint top-variant (the SNP with the lowest min-p across both traits), or on the GWAS lead consistently. Report PP under multiple centring choices; flag the locus if PP swings > 0.2 across centrings.

**Operational steps to diagnose window-centring bias** -- re-run coloc.abf three times with different window centres:

1. **GWAS-centred window:** +/- 500 kb around the GWAS lead SNP.
2. **eQTL-centred window:** +/- 500 kb around the eQTL top SNP for the gene.
3. **Joint top-variant window:** +/- 500 kb around the SNP with the lowest min-p across both traits.

Report all three PP.H4 values. If they agree within 0.1, the result is stable. If they swing > 0.2, the locus is borderline and the report must list all three centrings. For multi-causal loci (allelic heterogeneity), the joint top-variant window typically gives the most-defensible result for coloc.abf; coloc.susie removes the centring sensitivity by construction.
