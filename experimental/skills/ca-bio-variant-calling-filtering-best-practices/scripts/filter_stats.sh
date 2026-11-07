#!/bin/bash
# Compare a VCF before and after filtering: record counts, Ti/Tv, and counts per FILTER label.
#
# Usage:  filter_stats.sh <before.vcf[.gz]> <after.vcf[.gz]> [wgs|wes]
# Needs:  bcftools. Checked on bcftools 1.24.
# Warns:  Ti/Tv fell after filtering (the filters removed true transitions preferentially), or is below the
#         assay's floor (about 2.0 WGS, 3.0 WES: too loose), or the filter removed nothing.

set -euo pipefail

BEFORE=${1:-}
AFTER=${2:-}
ASSAY=${3:-}
if [ -z "$BEFORE" ] || [ -z "$AFTER" ]; then
    echo "Usage: filter_stats.sh <before.vcf[.gz]> <after.vcf[.gz]> [wgs|wes]" >&2
    exit 1
fi

n_before=$(bcftools view -H "$BEFORE" | wc -l)
n_after=$(bcftools view -H "$AFTER" | wc -l)
t_before=$(bcftools stats "$BEFORE" | grep '^TSTV' | cut -f5)
t_after=$(bcftools stats "$AFTER" | grep '^TSTV' | cut -f5)
echo "records: before $n_before, after $n_after, removed $((n_before - n_after))"
echo "Ti/Tv:   before ${t_before:-NA}, after ${t_after:-NA}  (expected about 2.0-2.1 WGS, 3.0-3.3 WES)"
echo "FILTER labels in $AFTER:"
bcftools query -f '%FILTER\n' "$AFTER" | sort | uniq -c

if [ "$n_after" -eq "$n_before" ]; then
    echo "NOTE: nothing was removed; the filter had no effect on these records." >&2
fi
if [ -n "$t_before" ] && [ -n "$t_after" ] && awk -v b="$t_before" -v a="$t_after" 'BEGIN{exit !(a < b)}'; then
    echo "WARNING: Ti/Tv fell after filtering ($t_before to $t_after): the filters are removing true transitions; relax them." >&2
fi
case "$ASSAY" in
    wgs) floor=2.0 ;;
    wes) floor=3.0 ;;
    *) floor= ;;
esac
if [ -n "$floor" ] && [ -n "$t_after" ] && awk -v t="$t_after" -v f="$floor" 'BEGIN{exit !(t < f)}'; then
    echo "WARNING: Ti/Tv $t_after is below $floor for $ASSAY: filtering is too loose." >&2
fi
