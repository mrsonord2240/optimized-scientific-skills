# Phenome-Wide Drug-Target MR

Moved verbatim from SKILL.md; the code block is now `scripts/phewas_curated_endpoints.R`. Read when scanning one target across many outcomes for on-target adverse effects.

## Phenome-Wide Drug-Target MR

**Goal:** For a single drug target (single protein), test causal effect across hundreds of outcomes to discover on-target adverse effects.

**Approach:** Hold cis-pQTL instrument set fixed; loop outcome over OpenGWAS catalogue or FinnGen DF12; multi-test correct over outcomes.

```bash
Rscript scripts/phewas_curated_endpoints.R pcsk9_cis_pqtls.tsv finngen_DF12_endpoints.tsv pcsk9_phewas_mr.tsv
```

The script takes the cis-pQTL table (SNP, BETA, SE, A1, A2, EAF, P), a curated endpoint table with an `id` column, an output path, and optional minimum sample size (default 50000) and population (default European). It needs an OpenGWAS token. `examples/phewas_drug_target_mr.R` is the fuller worked example (clumping, FDR, forest-plot input).

Document the curated endpoint list (FinnGen DF12, Open Targets curated trait map, or a manuscript-specific phecode hierarchy) in methods. The `outcomes_filt$id[1:200]` pattern is a debug shortcut, not a defensible pheWAS protocol. The PCSK9 -> T2D signal (Schmidt 2017 Lancet Diabetes Endocrinol 5:97) was discovered exactly via curated-endpoint pheWAS: cis-MR of LDL-lowering instruments revealed on-target T2D risk before clinical trials confirmed it. Drug-target pheWAS is the canonical use case.

