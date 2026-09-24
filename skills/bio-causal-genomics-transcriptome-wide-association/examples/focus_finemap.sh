#!/usr/bin/env bash
# Reference: MetaXcan 0.7+, FUSION 2.0+, FOCUS (pyfocus) 0.802, plink2 | Verify API if version differs
# FOCUS probabilistic gene-level fine-mapping of TWAS hits.
# Resolves co-significant gene clusters at gene-dense loci into a credible causal-gene set
# with per-gene posterior inclusion probabilities (PIPs).
#
# Install (checked 2026-09-21; a bare `pip install pyfocus` is non-functional -- see
# SKILL.md Tool Install Notes for the required post-install patches, including the two
# finemap.py fixes without which the run crashes at "Calculating PIPs"):
#   pip install pyfocus "pandas<2.2" "setuptools<81"
#
# Windows: use paths relative to the working directory below, never an absolute `F:/...`
# path -- pyfocus splits every positional argument on ':' to detect multi-ancestry input,
# which also splits a drive-letter colon (confirmed: mis-detects "2 populations" from 1 file).

set -euo pipefail

# ---- Inputs ----
GWAS_FILE='gwas.sumstats'                                # GWAS sumstats (CHR SNP BP A1 A2 Z P columns at minimum)
LD_REF_PREFIX='1000G_EUR/all'                            # One PLINK bfile (all.bim/bed/fam) covering every
                                                          # chromosome to analyze -- NOT chr-templated like
                                                          # FUSION's --ref_ld_chr. `focus finemap` passes this
                                                          # straight to pandas_plink.read_plink(); a literal
                                                          # ".../chr" prefix with no matching file 404s
                                                          # (confirmed 2026-09-21). --chr/--locations subset
                                                          # the loaded SNPs afterward, they do not pick a file.
                                                          # For real per-chromosome files, use pandas_plink's
                                                          # own glob syntax instead, e.g. '1000G_EUR/chr*.bed'.
FOCUS_DB='focus_gtex_v8_whole_blood.db'                  # FOCUS DB matched to the TWAS weight panel
TISSUE='Whole_Blood'                                     # Tissue label inside the FOCUS DB
LOCATIONS='38:EUR'                                       # Required (not optional) in installed pyfocus 0.802;
                                                          # use 37:EUR instead for GRCh37-aligned panels

# Genome-wide significance threshold for SNPs that flag a locus for fine-mapping.
# 5e-8 is standard GWAS genome-wide significance; loci below this threshold trigger FOCUS.
P_THRESHOLD='5e-8'

OUT_PREFIX='gwas_focus_whole_blood'

# ---- Step 1: build FOCUS database from FUSION weights if not pre-built ----
# Skip this if using a pre-built FOCUS DB. Custom panels need this step; `focus import` needs
# mygene + rpy2 (see SKILL.md, FOCUS section, for the direct-build alternative).
# focus import gtex_whole_blood.pos fusion --tissue Whole_Blood --output focus_gtex_v8_whole_blood

# ---- Step 2: run FOCUS fine-mapping ----
# FOCUS takes the underlying GWAS sumstats (not the TWAS Z output) and internally:
#   1. computes per-gene TWAS Z using the panel weights
#   2. estimates the gene-by-gene predicted-expression correlation matrix as the LD analog
#   3. runs a Bayesian variable-selection model over genes (analog of variant fine-mapping)
focus finemap \
    "${GWAS_FILE}" \
    "${LD_REF_PREFIX}" \
    "${FOCUS_DB}" \
    --p-threshold "${P_THRESHOLD}" \
    --tissue "${TISSUE}" \
    --locations "${LOCATIONS}" \
    --out "${OUT_PREFIX}"
# --locations is required here, not optional: omitting it crashes with "Please specify
# independent regions location or default regions with '37:EUR', etc." (confirmed
# 2026-09-19; there is no default-LD-block fallback in installed pyfocus 0.802). The same
# build:pop syntax extends to multi-ancestry 38:EUR-EAS-AFR (with colon-separated
# per-ancestry sumstats/LD/weight DBs) for MA-FOCUS (mancusolab/ma-focus) below.

# ---- Step 3: filter credible-set genes ----
# FOCUS reports per-gene PIP. PIP >= 0.8 = causal candidate; 0.5 <= PIP < 0.8 = suggestive;
# PIP < 0.5 at a co-significant locus = LD-tagged co-regulated gene (NOT causal).
# Mancuso 2019 Nat Genet 51:675 convention.
#
# Installed pyfocus 0.802 never writes a plain "pip" column: single-ancestry output names it
# "pips_pop1" (population-indexed even for one population); multi-ancestry adds "pips_me" for
# the cross-ancestry marginal PIP (confirmed from pyfocus/finemap.py's create_output(), 2026-09-19).
PIP_THRESH='0.8'
PIP_COL='pips_pop1'   # use 'pips_me' for the cross-ancestry PIP on MA-FOCUS output

awk -v t="${PIP_THRESH}" -v pcol="${PIP_COL}" -F'\t' '
    NR==1 {print; for (i=1; i<=NF; i++) col[$i]=i; next}
    $col[pcol] >= t {print}
' "${OUT_PREFIX}.focus.tsv" > "${OUT_PREFIX}.credible.tsv"

echo 'FOCUS fine-mapping complete.'
echo "Full output: ${OUT_PREFIX}.focus.tsv"
echo "Causal-gene candidates (PIP >= ${PIP_THRESH}): ${OUT_PREFIX}.credible.tsv"

# ---- Step 4 (optional): MA-FOCUS for multi-ancestry ----
# Requires per-ancestry sumstats, per-ancestry LD references, and per-ancestry FOCUS DBs.
# Joint inference assumes a shared causal gene across ancestries; trans-ethnic
# gene-effect heterogeneity violates this and produces inflated heterogeneous-group probability.
#
# MA-FOCUS reuses the single-ancestry `focus finemap` CLI; multi-ancestry mode
# is signaled by colon-separated per-ancestry sumstats/LD/weight DBs plus
# paired ancestry codes in --locations.
# focus finemap \
#     eur.sumstats.tsv.gz:eas.sumstats.tsv.gz:afr.sumstats.tsv.gz \
#     1000G_EUR/chr:1000G_EAS/chr:1000G_AFR/chr \
#     focus_eur.db:focus_eas.db:focus_afr.db \
#     --p-threshold "${P_THRESHOLD}" \
#     --tissue "${TISSUE}" \
#     --locations 38:EUR-EAS-AFR \
#     --out gwas_ma_focus
