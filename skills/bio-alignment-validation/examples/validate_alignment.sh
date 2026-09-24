#!/bin/bash
# Reference: samtools 1.19+ (checked on samtools 1.24) | Verify API if version differs
# Comprehensive alignment validation with a pass/warn/fail verdict.
#
# Usage: validate_alignment.sh <bam_file>          (index needed only for the idxstats section)
# Exit status: 0 = every metric PASS or WARN, 1 = at least one metric in the FAIL band,
#              2 = file missing / empty / unreadable / fails samtools quickcheck / no primary records.
# Every rate counts primary records only (-F 2304; mapped: -F 2308), unplaced unmapped reads included.
# Bands: see SKILL.md "Quality Thresholds Summary" (germline short-read DNA defaults).

set -uo pipefail

BAM=${1:-}
if [ -z "$BAM" ]; then
    echo "Usage: $0 <bam_file>" >&2
    exit 2
fi
if ! { test -s "$BAM" && samtools quickcheck -v "$BAM"; }; then
    echo "FAIL: $BAM is missing, empty or fails samtools quickcheck" >&2
    exit 2
fi

NAME=$(basename "$BAM" .bam)
FAILS=""
WARNS=""

# record <label> <PASS|WARN|FAIL>: print the grade and remember WARN/FAIL for the verdict
record() {
    echo "  -> $1: $2"
    [ "$2" = FAIL ] && FAILS="${FAILS:+$FAILS, }$1"
    [ "$2" = WARN ] && WARNS="${WARNS:+$WARNS, }$1"
    return 0
}
# grade <label> <value> <good_above> <warn_from>  (higher is better)
grade() {
    record "$1" "$(awk -v v="$2" -v good="$3" -v warn="$4" 'BEGIN{print (v>good)?"PASS":(v>=warn)?"WARN":"FAIL"}')"
}

# Decode primary records once. Unlike `view -c`, this also fails when a CRAM's
# reference cannot be resolved; keep all counters in this one awk pass.
if ! metrics=$(samtools view -F 2304 "$BAM" | awk '
function hasbit(value, bit) { return int(value / bit) % 2 }
{
    total++
    flag = $2
    if (hasbit(flag, 4)) next
    mapped++
    if (hasbit(flag, 1)) paired++
    if (hasbit(flag, 2)) proper++
    if (hasbit(flag, 16)) reverse++; else forward++
    mapq_sum += $5; mapq_count++
    if ($5 >= 30) mapq_ge_30++
    if ($5 >= 1) mapq_ge_1++
}
END { printf "%d\t%d\t%d\t%d\t%d\t%d\t%.12f\t%d\t%d\t%d\n", total, mapped, paired, proper, forward, reverse, mapq_sum, mapq_count, mapq_ge_30, mapq_ge_1 }
'); then
    echo "FAIL: $BAM cannot be decoded completely (supply REF_PATH or -T for CRAM)" >&2
    exit 2
fi
IFS=$'\t' read -r total mapped paired proper fwd rev mapq_sum mapq_count q30 q1 <<< "$metrics"
if [ "$total" -eq 0 ]; then
    echo "FAIL: $BAM has no primary records; nothing to validate" >&2
    exit 2
fi

echo "=== Alignment Validation: $NAME ==="

echo -e "\n--- Basic Stats ---"
samtools flagstat "$BAM"

echo -e "\n--- Mapping Rate ---"
map_rate=$(awk -v a="$mapped" -v b="$total" 'BEGIN{printf "%.2f", 100*a/b}')
echo "Mapped: $mapped / $total primary records (${map_rate}%)"
grade "Mapping rate" "$map_rate" 95 90

echo -e "\n--- Proper Pairing ---"
if [ "$total" -lt 100 ]; then
    echo "Too few primary reads to grade pairing or strand balance (n<100)"
elif [ "$paired" -gt 0 ]; then
    pair_rate=$(awk -v a="$proper" -v b="$paired" 'BEGIN{printf "%.2f", 100*a/b}')
    echo "Properly paired: $proper / $paired mapped paired reads (${pair_rate}%)"
    grade "Proper pairing" "$pair_rate" 90 80
else
    echo "No mapped paired reads (single-end or unpaired): proper pairing not graded"
fi

echo -e "\n--- Insert Size ---"
samtools stats "$BAM" 2>/dev/null | grep '^SN' | grep 'insert size' | cut -f2,3

if [ "$mapped" -gt 0 ] && [ "$total" -ge 100 ]; then
    echo -e "\n--- Strand Balance ---"
    fwd=$(samtools view -c -F 2324 "$BAM")
    rev=$(samtools view -c -f 16 -F 2308 "$BAM")
    strand=$(awk -v f="$fwd" -v r="$rev" 'BEGIN{printf "%.3f", f/(f+r)}')
    echo "Forward: $fwd, Reverse: $rev, Forward fraction F/(F+R): $strand"
    record "Strand balance" "$(awk -v s="$strand" 'BEGIN{print (s>=0.48&&s<=0.52)?"PASS":(s>=0.45&&s<=0.55)?"WARN":"FAIL"}')"

    echo -e "\n--- MAPQ (mapped primary reads) ---"
    echo "MAPQ 0 (multi-mapper): $((mapped - q1))"
    echo "MAPQ >= 30: $q30 ($(awk -v a="$q30" -v b="$mapped" 'BEGIN{printf "%.1f", 100*a/b}')%)"
    mean_mapq=$(awk -v s="$mapq_sum" -v n="$mapq_count" 'BEGIN{printf "%.1f", s/n}')
    echo "Mean MAPQ: $mean_mapq"
    grade "Mean MAPQ" "$mean_mapq" 40 30
fi

echo -e "\n--- Top Chromosomes ---"
if idx=$(samtools idxstats "$BAM" 2>/dev/null); then
    echo "$idx" | awk '$1!="*"' | sort -k3,3nr | head -10
else
    echo "(no index: run samtools index to list per-chromosome counts)"
fi

echo -e "\n--- Verdict ---"
[ -n "$FAILS" ] && echo "FAIL: $FAILS"
[ -n "$WARNS" ] && echo "WARN: $WARNS"
if [ -n "$FAILS" ]; then
    exit 1
fi
[ -z "$WARNS" ] && echo "All metrics within normal range"
exit 0
