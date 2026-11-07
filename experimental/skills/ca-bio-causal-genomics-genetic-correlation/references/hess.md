# rho-HESS: Alternative Local rg
**Goal:** Per-locus bivariate rg using the HESS quadratic-form estimator; complementary to LAVA.

**Approach:** Estimate local h2 per trait first; then bivariate cross-trait estimator using LD-projected effect estimates per locus.

```bash
# Step 1: eigenvalues + projections for both traits (per chromosome; two sumstats SPACE-separated).
# --local-rhog auto-writes per-trait files (step1_trait1_*, step1_trait2_*) AND the covariance
# intermediates (step1_chrN.eig.gz, step1_chrN.prjprod.gz) under the shared --out prefix.
for chr in {1..22}; do
    hess.py \
        --local-rhog trait1.sumstats.gz trait2.sumstats.gz \
        --chrom $chr \
        --bfile 1kg_EUR_chr${chr} \
        --partition fourier_ls-chr${chr}.bed \
        --out step1
done

# Step 2: per-trait local h2 from the Step-1 outputs (MUST run before Step 3)
hess.py --prefix step1_trait1 --out step2_trait1
hess.py --prefix step1_trait2 --out step2_trait2

# Step 3: local genetic covariance / rg. --local-hsqg-est passes the Step-2 per-trait local h2,
# --num-shared is the overlap count, --pheno-cor the phenotypic correlation (any value when
# --num-shared 0). Auto-aggregates across all chromosomes; no per-chromosome loop here.
hess.py \
    --prefix step1 \
    --local-hsqg-est step2_trait1.txt step2_trait2.txt \
    --num-shared 0 \
    --pheno-cor 0 \
    --out step3
```

`--num-shared` is the number of overlapping individuals; set 0 only if truly independent. HESS partition files use the Berisa & Pickrell 2016 LD-block boundaries (`fourier_ls-*.bed`). LAVA has largely superseded HESS for new analyses, but rho-HESS remains in active use for replication / triangulation.
