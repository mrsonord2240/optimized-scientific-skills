#!/usr/bin/env bash
# Reference: MAGMA 1.10+, PLINK 1.9+, 1000 Genomes EUR reference | Verify CLI flags if version differs
#
# MAGMA gene-based and gene-set analysis from GWAS summary statistics.
# Three stages:
#   1. SNP-to-gene annotation (window-based; FUMA default 35kb upstream + 10kb downstream)
#   2. Gene-based association (multi-SNP joint test per gene; LD-corrected via reference panel)
#   3. Gene-set enrichment (competitive test against MSigDB C2)

set -euo pipefail

GWAS_PVAL='gwas.pval.tsv'
GWAS_NCOL='N'
GENE_LOC='NCBI37.3.gene.loc'
SNP_LOC='gwas.snploc'
REF_BFILE='g1000_eur'
GENESET_GMT='msigdb_v7.5_C2.gmt'
OUT_PREFIX='magma_run'

# Window: 35 kb upstream + 10 kb downstream (FUMA recommendation; balances regulatory capture vs gene-dense dilution)
WINDOW='35,10'

# Step 1: SNP-to-gene annotation
magma --annotate window="$WINDOW" \
    --snp-loc "$SNP_LOC" \
    --gene-loc "$GENE_LOC" \
    --out "${OUT_PREFIX}_annot"

# Step 2: Gene-based association
magma --bfile "$REF_BFILE" \
    --pval "$GWAS_PVAL" "ncol=${GWAS_NCOL}" \
    --gene-annot "${OUT_PREFIX}_annot.genes.annot" \
    --out "${OUT_PREFIX}_gene"

# Windows note: MAGMA 1.10 on Windows writes "<prefix>.genes.out.txt" (extra .txt)
# instead of "<prefix>.genes.out" -- the same mismatch documented in
# references/pops.md for pops.py. Resolve it here too, since the
# Bonferroni/top-50 steps below read .genes.out directly (verified: this script fails
# with "No such file or directory" on Windows without this).
GENES_OUT="${OUT_PREFIX}_gene.genes.out"
if [ ! -f "$GENES_OUT" ] && [ -f "${GENES_OUT}.txt" ]; then
    cp "${GENES_OUT}.txt" "$GENES_OUT"
fi

# Step 3: Gene-set enrichment (competitive test; preferred over self-contained)
# MINIMUM GENE COUNT: this regression conditions on 6 internal covariates (gene size,
# log(gene size), gene density, log(gene density), inverse MAC, log(inverse MAC)).
# With too few genes in --gene-results relative to those 6 covariates, the regression
# is structurally unidentifiable and MAGMA aborts with "ERROR: insufficient degrees of
# freedom to run analyses" (or, with even fewer genes, "no input variables to analyse"
# from a gene-set variance check that runs first) -- verified on real locus-scale runs
# with 3 and 5 genes. Rule of thumb: gene-set enrichment needs several hundred+ genes;
# a single-locus run (a handful of genes) should SKIP this step entirely -- gene-set
# enrichment is a genome-wide- or many-loci-scale analysis, not a per-locus one.
N_GENES_STEP3=$(grep -vc '^#' "${OUT_PREFIX}_gene.genes.raw")
if [ "$N_GENES_STEP3" -lt 200 ]; then
    echo "Skipping Step 3 (gene-set enrichment): only ${N_GENES_STEP3} genes in" \
         "${OUT_PREFIX}_gene.genes.raw, well under the several-hundred+ this competitive" \
         "regression needs. Run gene-set enrichment only on genome-wide (or many-loci) MAGMA output."
else
    magma --gene-results "${OUT_PREFIX}_gene.genes.raw" \
        --set-annot "$GENESET_GMT" \
        --out "${OUT_PREFIX}_geneset"
fi

# Outputs:
#   ${OUT_PREFIX}_gene.genes.out: per-gene Z, p, n_SNPs
#   ${OUT_PREFIX}_geneset.gsa.out: per-set p, beta, beta_se

# Bonferroni threshold across the genes MAGMA tested. .genes.out has a header row, so skip
# it when counting. awk does the division (portable: `bc` is absent from Git-for-Windows bash).
N_GENES=$(tail -n +2 "$GENES_OUT" | wc -l | tr -d ' ')
echo "Bonferroni threshold at alpha=0.05 across ${N_GENES} genes: $(awk -v n="$N_GENES" 'BEGIN{printf "%.6g", 0.05/n}')"

# Top 50 genes by p
# MAGMA .genes.out columns: GENE CHR START STOP NSNPS NPARAM N ZSTAT P
# P is the 9th whitespace-separated field; sort ascending numerically.
{ head -1 "$GENES_OUT"
  tail -n +2 "$GENES_OUT" | sort -k9,9g
} | head -51 > "${OUT_PREFIX}_top50.tsv"

echo "MAGMA gene-based and gene-set complete."
echo "Pair with PoPS (FinucaneLab/pops): pass --magma_prefix ${OUT_PREFIX}_gene (the prefix, not a file name)."
