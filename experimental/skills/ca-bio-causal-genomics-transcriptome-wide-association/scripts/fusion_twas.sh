#!/usr/bin/env bash
# fusion_twas.sh -- FUSION sumstat TWAS per chromosome, then conditional/joint analysis.
#
# Inputs (positional):
#   SUMSTATS     GWAS sumstats with columns SNP A1 A2 Z (Z on the standardised scale)
#   WEIGHTS_POS  FUSION .pos file for the weight panel (gusevlab.org/projects/fusion)
#   WEIGHTS_DIR  directory holding the per-gene .wgt.RDat files named in the .pos file
#   LD_PREFIX    ancestry-matched LD reference prefix, e.g. 1000G_EUR_LD/EUR.  (files EUR.<chr>.bed/.bim/.fam)
#   OUT_PREFIX   output prefix: <prefix>_chr<N>.dat, <prefix>_all.dat, <prefix>_joint_chr<N>.*
#   CHRS         optional, space-separated chromosome list (default: "1 2 ... 22")
# Environment:
#   RSCRIPT      R launcher (default: Rscript). Point it at a wrapper if bare Rscript is not usable.
#   LOCUS_WIN    conditioning window in bp for post_process (default 100000; conservative for non-HLA loci)
# Run from inside the fusion_twas clone (FUSION sources utils/ through here(), so its scripts must be
# run from that directory); give the input and output paths relative to it, or absolute.
# Usage (cwd = fusion_twas clone; the OUT_PREFIX directory must exist):
#   /path/to/scripts/fusion_twas.sh gwas.sumstats gtex_whole_blood.pos gtex_whole_blood_wgt/ \
#       1000G_EUR_LD/EUR. twas_out/twas "$(seq -s' ' 1 22)"
# post_process.R runs once per chromosome that has TWAS output; twas_joint*.dat carries per-gene
# conditional Z (not PIPs; use FOCUS for those). On single-SNP ("top1") weights, post_process.R needs
# the drop=FALSE patch described in SKILL.md Common Errors.
set -euo pipefail

SUMSTATS=${1:?SUMSTATS}; WEIGHTS_POS=${2:?WEIGHTS_POS}; WEIGHTS_DIR=${3:?WEIGHTS_DIR}
LD_PREFIX=${4:?LD_PREFIX}; OUT_PREFIX=${5:?OUT_PREFIX}
CHRS=${6:-$(seq -s' ' 1 22)}
RSCRIPT=${RSCRIPT:-Rscript}; LOCUS_WIN=${LOCUS_WIN:-100000}

for chr in ${CHRS}; do
    "${RSCRIPT}" FUSION.assoc_test.R \
        --sumstats "${SUMSTATS}" \
        --weights "${WEIGHTS_POS}" \
        --weights_dir "${WEIGHTS_DIR}" \
        --ref_ld_chr "${LD_PREFIX}" \
        --chr "${chr}" \
        --out "${OUT_PREFIX}_chr${chr}.dat"
done
# Keep one header line when concatenating the per-chromosome tables
first=1
: > "${OUT_PREFIX}_all.dat"
for chr in ${CHRS}; do
    f="${OUT_PREFIX}_chr${chr}.dat"
    [ -s "${f}" ] || continue
    if [ "${first}" = 1 ]; then cat "${f}" >> "${OUT_PREFIX}_all.dat"; first=0; else tail -n +2 "${f}" >> "${OUT_PREFIX}_all.dat"; fi
done

# Conditional joint analysis per chromosome; joint Z > 4 marks an independent gene
for chr in ${CHRS}; do
    [ -s "${OUT_PREFIX}_chr${chr}.dat" ] || continue
    "${RSCRIPT}" FUSION.post_process.R \
        --sumstats "${SUMSTATS}" \
        --input "${OUT_PREFIX}_all.dat" \
        --out "${OUT_PREFIX}_joint_chr${chr}.dat" \
        --ref_ld_chr "${LD_PREFIX}" \
        --chr "${chr}" \
        --plot --locus_win "${LOCUS_WIN}"
done
