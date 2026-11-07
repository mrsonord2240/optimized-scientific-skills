#!/bin/bash
# Purpose: normalize a VCF, then annotate with Ensembl VEP on MANE Select with a constrained pick order.
# Usage:   scripts/annotate_vep.sh input.vcf.gz reference.fa out_prefix [--cache-dir DIR] [--assembly GRCh38]
#            [--fork 4] [--plugin "dbNSFP,dbNSFP4.3a.gz,REVEL_score" ...]
#          --plugin may repeat; give ONE calibrated predictor per evidence type (see routes/predictors.md).
# Output:  out_prefix_norm.vcf.gz, out_prefix_vep.vcf.gz (+ index), out_prefix_vep.summary.txt
#          (VEP version, assembly, record count, records with a MANE Select block, HIGH/MODERATE count).
# Refuses: an assembly other than GRCh38 (MANE exists only there), a missing input, no vep on PATH, an output
#          whose CSQ header has no MANE_SELECT field.
set -euo pipefail

die() { echo "REFUSED: $*" >&2; exit 2; }
[ $# -ge 3 ] || die "usage: $0 input.vcf.gz reference.fa out_prefix [--cache-dir DIR] [--assembly GRCh38] [--fork N] [--plugin SPEC]..."
INPUT=$1; REFERENCE=$2; OUT=$3; shift 3
CACHE=""; ASSEMBLY="GRCh38"; FORK=4; PLUGINS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --cache-dir) CACHE=$2; shift 2;;
    --assembly) ASSEMBLY=$2; shift 2;;
    --fork) FORK=$2; shift 2;;
    --plugin) PLUGINS+=(--plugin "$2"); shift 2;;
    *) die "unknown option $1";;
  esac
done
[ "$ASSEMBLY" = "GRCh38" ] || die "MANE Select exists only on GRCh38; for $ASSEMBLY lift over or pin a per-gene transcript list and annotate by hand"
command -v bcftools >/dev/null 2>&1 || die "bcftools not found on PATH"
command -v vep >/dev/null 2>&1 || die "vep not found on PATH"
[ -f "$INPUT" ] || die "file not found: $INPUT"
[ -f "$REFERENCE" ] || die "reference not found: $REFERENCE"
bcftools index -s "$INPUT" >/dev/null 2>&1 || die "indexed VCF required: $INPUT"

bcftools norm -f "$REFERENCE" -m-any "$INPUT" -Oz -o "${OUT}_norm.vcf.gz"
bcftools index -f "${OUT}_norm.vcf.gz"

CACHE_ARGS=()
[ -z "$CACHE" ] || CACHE_ARGS=(--dir_cache "$CACHE")
# --mane_select with a pick order where MANE leads, so transcript length or accession never decides the reported block.
vep -i "${OUT}_norm.vcf.gz" -o "${OUT}_vep.vcf" --vcf --force_overwrite --cache --offline ${CACHE_ARGS[@]+"${CACHE_ARGS[@]}"} \
    --assembly "$ASSEMBLY" --everything --mane_select --pick \
    --pick_order mane_select,mane_plus_clinical,canonical,biotype,rank --fork "$FORK" ${PLUGINS[@]+"${PLUGINS[@]}"}

grep -q '^##INFO=<ID=CSQ,.*MANE_SELECT' "${OUT}_vep.vcf" || die "VEP output has no MANE_SELECT field; the cache may predate MANE (VEP 110+ needed)"
bgzip -f "${OUT}_vep.vcf"
bcftools index -f "${OUT}_vep.vcf.gz"

{
  echo "vep: $(vep --help 2>&1 | grep -m1 -i 'ensembl-vep' || echo unknown)"
  echo "assembly: $ASSEMBLY"
  echo "records: $(bcftools view -H "${OUT}_vep.vcf.gz" | wc -l | tr -d ' ')"
  echo "records with a MANE Select block: $(bcftools view -H "${OUT}_vep.vcf.gz" | grep -c 'MANE_SELECT\|NM_[0-9]*\.[0-9]*' || true) (approximate: counts INFO text carrying a MANE accession)"
  echo "records with HIGH or MODERATE impact (triage, not classification): $(bcftools view -H "${OUT}_vep.vcf.gz" | grep -c '|HIGH|\||MODERATE|' || true)"
} > "${OUT}_vep.summary.txt"
cat "${OUT}_vep.summary.txt"
echo "wrote ${OUT}_vep.vcf.gz and ${OUT}_vep.summary.txt"
