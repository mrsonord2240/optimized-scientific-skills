#!/usr/bin/env bash
# Reference: DAP-G (xqwen/dap, built from source), PLINK 1.9+ | Verify CLI flags with `dap-g -d_z --help`-style docs (repo README) if version differs
##
## DAP-G pipeline for a single GWAS locus, sharing the same locus-extraction
## and in-sample-LD steps as examples/finemap_pipeline.sh -- only the input
## file shape and the run/report step differ.
## Requires: dap-g binary, plink, GWAS summary stats, reference panel.

set -euo pipefail

# --- Configuration ---
GWAS_FILE="gwas_sumstats.txt"    # Columns: SNP CHR POS A1 A2 MAF BETA SE
REF_PANEL="1000G_EUR"            # Plink binary prefix for LD reference
CHR=6
START=30000000
END=31000000

# --- Step 1: Extract locus SNPs (same as examples/finemap_pipeline.sh) ---
awk -v chr="$CHR" -v start="$START" -v end="$END" \
  '$2 == chr && $3 >= start && $3 <= end' "$GWAS_FILE" > locus_gwas.txt
awk '{print $1}' locus_gwas.txt > locus_snps.txt

# --- Step 2: In-sample LD, space-delimited (DAP-G, like FINEMAP, needs a
#     plain space-separated matrix; PLINK's --r square default is tab-delimited) ---
plink --bfile "$REF_PANEL" \
  --chr "$CHR" --from-bp "$START" --to-bp "$END" \
  --extract locus_snps.txt \
  --r square spaces \
  --out locus_ld
mv locus_ld.ld locus.LD.dat

# --- Step 3: Build DAP-G's z-value file (snp_id z_value, no header) ---
# Columns in locus_gwas.txt: SNP CHR POS A1 A2 MAF BETA SE
awk '{print $1, $7/$8}' locus_gwas.txt > locus.zval.dat

# --- Step 4: Run DAP-G ---
# Exits with code 1 even on success -- judge by stdout, not the exit code
# (confirmed on this build: the package's own bundled example also exits 1).
set +e
dap-g -d_z locus.zval.dat -d_ld locus.LD.dat -t 4 > locus.dap.out 2>&1
set -e

# --- Step 5: Report results ---
echo ""
echo "=== DAP-G Independent Association Signal Clusters ==="
grep -A 20 "Independent association signal clusters" locus.dap.out || echo "(no clusters found; check locus.dap.out)"

echo ""
echo "Top 10 SNPs by PIP:"
grep -E '^\(\(' locus.dap.out | sort -t$'\t' -k3 -gr | head -10
