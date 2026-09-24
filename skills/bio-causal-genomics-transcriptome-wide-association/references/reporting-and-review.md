# Required Reporting and Reviewer Pushback

## Required Reporting for Publication

A defensible TWAS report includes every item below in methods or supplement:

- GWAS sumstat source, effective N (Neff), and ancestry composition
- Weight panel and version (e.g. "GTEx v8 MASHR-EUR, PredictDB release 2022-01")
- LD reference panel and version (e.g. "1000 Genomes Phase 3 EUR" or "UK Biobank array")
- Per-tissue list (or explicit "all 49 GTEx v8 tissues")
- Multiple-testing correction strategy (S-MultiXcan joint, per-tissue Bonferroni, or per-tissue FDR)
- FOCUS PIP threshold and prior-probability sensitivity scan
- Coloc PP.H4 threshold (cross-reference causal-genomics/colocalization-analysis)
- cis-eQTL MR estimate, instrument F-statistic, and TwoSampleMR package version (cross-reference causal-genomics/mendelian-randomization)
- Triangulation rule (e.g. "3-of-4 concordance across TWAS, FOCUS, coloc, cis-MR")
- HLA exclusion confirmed (chr6:25-35 Mb dropped from genome-wide summaries)

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "LD-induced false positive at gene-dense locus?" | FOCUS PIP reported per gene; only PIP >= 0.8 carried forward as candidate causal |
| "Tissue mis-specification?" | Tissue Selection Protocol applied (`tissue-and-model-selection.md`) (sLDSC + CELLEX + MAGMA); cross-tissue S-MultiXcan reported for replication |
| "Why GTEx and not eQTLGen?" | Justified by trait biology: blood-relevant traits use eQTLGen (N ~ 31k); tissue-specific traits use GTEx v8 MASHR |
| "MHC?" | chr6:25-35 Mb excluded; HLA-TAPAS run separately for HLA-relevant traits |
| "Triangulation?" | TWAS + coloc + cis-MR + FOCUS run; 3-of-4 concordance required for the strong-candidate label |
| "Ancestry transfer?" | Ancestry-matched weights used where available; MA-FOCUS applied for multi-ancestry GWAS |
| "Prior sensitivity in FOCUS?" | `--prior-prob` scanned at 1e-2, 1e-3, 1e-4; PIPs reported across the scan |
