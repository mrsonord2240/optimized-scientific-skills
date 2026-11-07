#!/bin/bash
# Extra thresholds on a FilterMutectCalls VCF, tested on the TUMOR column only.
#
# Usage:  somatic_filter.sh <mutect2_filtered.vcf[.gz]> <out.vcf> [tlod=6.3] [af=0.05] [dp=20] [tumor=<sample>] [pass_only=0]
#           tumor=       tumor sample name; default is the ##tumor_sample header line
#           pass_only=1  also require FILTER=PASS (drop what FilterMutectCalls rejected). Default 0 keeps the
#                        source behaviour: the thresholds are applied to every record whatever its FILTER.
# Needs:  bcftools. Checked on bcftools 1.24.
# Traps:  exit 2 when the tumor sample cannot be found (without this guard bcftools segfaults, exit 139, on a
#         tumor index of -1). Mutect2 often lists the normal first, so FMT/AF[0] would test the normal.

set -euo pipefail

IN=${1:-}
OUT=${2:-}
if [ -z "$IN" ] || [ -z "$OUT" ]; then
    echo "Usage: somatic_filter.sh <in.vcf[.gz]> <out.vcf> [tlod=6.3] [af=0.05] [dp=20] [tumor=NAME] [pass_only=0]" >&2
    exit 1
fi
shift 2
TLOD=6.3; AF=0.05; DP=20; TUMOR=; PASS_ONLY=0
for kv in "$@"; do
    case "$kv" in
        tlod=*) TLOD=${kv#tlod=} ;;
        af=*) AF=${kv#af=} ;;
        dp=*) DP=${kv#dp=} ;;
        tumor=*) TUMOR=${kv#tumor=} ;;
        pass_only=*) PASS_ONLY=${kv#pass_only=} ;;
        *) echo "ERROR: unknown argument $kv" >&2; exit 1 ;;
    esac
done

if [ -z "$TUMOR" ]; then
    TUMOR=$(bcftools view -h "$IN" | grep '^##tumor_sample=' | cut -d= -f2 || true)
fi
T=-1
if [ -n "$TUMOR" ]; then
    idx=$(bcftools query -l "$IN" | grep -nxF "$TUMOR" | cut -d: -f1 || true)
    if [ -n "$idx" ]; then T=$((idx - 1)); fi
fi
if [ -z "$TUMOR" ] || [ "$T" -lt 0 ]; then
    echo "ERROR: tumor sample '${TUMOR}' not found among: $(bcftools query -l "$IN" | tr '\n' ' '). Pass tumor=<name>." >&2
    exit 2
fi

EXPR="INFO/TLOD>$TLOD && FMT/AF[$T:0]>$AF && FMT/DP[$T]>$DP"
if [ "$PASS_ONLY" = 1 ]; then EXPR="FILTER=\"PASS\" && $EXPR"; fi
echo "tumor sample: $TUMOR (column $T); expression: $EXPR"
bcftools filter -i "$EXPR" "$IN" -o "$OUT"
echo "records: in $(bcftools view -H "$IN" | wc -l), out $(bcftools view -H "$OUT" | wc -l)"
