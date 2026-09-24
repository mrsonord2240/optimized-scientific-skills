# coloc.susie: LD matrix and multi-causal pipeline

## LD Matrix Construction for coloc.susie

**Requirements:**

- Signed Pearson r (not r2). coloc.susie expects directional LD; squared LD silently inverts effect-direction inference.
- Ancestry-matched to GWAS / eQTL ancestry. EUR LD on AFR z-scores produces spurious credible sets.
- SNP-order-aligned to the beta vector and named to match (row/column names = SNP IDs).
- Positive semi-definite. Numerical-noise negative eigenvalues must be repaired.
- Effective N sample-size-matched to the trait being fine-mapped (provide via `runsusie(..., n = N)`).

```bash
# plink2 phased r (signed Pearson); square matrix output
plink2 --pfile 1KG_EUR \
    --extract snps.txt --chr 6 --from-bp X --to-bp Y \
    --r-phased square --out locus_ld
```

```r
# Alternative: in-sample LD from BED via bigsnpr
R <- bigsnpr::snp_cor(snp_obj$genotypes, ind.col = locus_snps)
# PSD repair if negative eigenvalues from numerical noise
R <- as.matrix(Matrix::nearPD(R)$mat)
dimnames(R) <- list(snp_ids, snp_ids)
```

**Critical:** Row and column order of R MUST match SNP order in the beta vector -- silent failure otherwise. The SuSiE objective stays finite under mis-ordering and returns nonsense credible sets. Verify with `stopifnot(rownames(R) == names(beta))` before `runsusie`. Cross-reference causal-genomics/fine-mapping for the full LD diagnostic block (`estimate_s_rss` lambda < 0.05, `kriging_rss` outlier inspection).

## coloc.susie Multi-Causal Pipeline

**Goal:** Test colocalization at a locus with multiple independent signals (allelic heterogeneity).

**Approach:** Run SuSiE on each trait's summary statistics with ancestry-matched LD; verify LD-z-score consistency; coloc-test each pair of credible sets.

```bash
Rscript scripts/coloc_susie.R gwas.tsv eqtl.tsv ld.tsv --gwas-type cc --gwas-s 0.3 --gwas-n 50000     --eqtl-type quant --eqtl-sdy 1 --eqtl-n 500 --L 10 --out coloc_susie_out
```

`scripts/coloc_susie.R` runs the MHC / chr 8 gate (`scripts/flag_excluded_region.R`), checks SNP order across the two tables and the LD matrix, independently checks GWAS and eQTL `estimate_s_rss` lambda values against `--lambda-max` (default 0.05), runs `runsusie` on each trait with the same LD, then `coloc.susie`; each row of the result is one (hit1, hit2) credible-set pair. `ld.tsv` is the signed r matrix with SNP ids as first column and header.

LD matrix MUST be in the same SNP order as the beta vector; mis-ordering silently produces nonsense.

### coloc.susie -- LD reference mismatch

**Trigger:** Z-scores from GWAS / eQTL of one ancestry, LD matrix from 1000 Genomes EUR (or any non-matched reference).

**Mechanism:** SuSiE assumes z-scores and the supplied LD are jointly consistent. Ancestry mismatch or sample-size mismatch produces a non-positive-definite implicit covariance; SuSiE responds by returning spurious credible sets that include LD-mismatched SNPs.

**Symptom:** `susieR::estimate_s_rss(z, R, n)` returns lambda > 0.05; `susieR::kriging_rss` flags off-diagonal SNPs with extreme studentized residuals; credible sets are oddly large (50+ SNPs) or include SNPs distant in LD from the lead.

**Fix:** Use in-sample LD when at all possible (per-cohort plink `--r square`). If reference must be external, match ancestry (1KG superpopulation) and superpopulation-stratify. Run `estimate_s_rss` and report lambda; if > 0.05, drop the locus or switch to coloc.abf.
