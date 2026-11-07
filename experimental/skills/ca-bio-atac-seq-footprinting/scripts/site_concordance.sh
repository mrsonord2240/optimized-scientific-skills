#!/bin/bash
# Site-level concordance between TOBIAS motif-site calls and a second call set (HINT-ATAC or Wellington footprints, ChIP peaks).
# Usage: site_concordance.sh bound.bed unbound.bed other_calls.bed [regions.bed]
#   bound.bed / unbound.bed   TOBIAS BINDetect beds (<outdir>/bindetect/<TF>/beds/<TF>_<cond>_{bound,unbound}.bed; cat several TFs if wanted)
#   other_calls.bed           footprints or peaks from the second method
#   regions.bed (optional)    restrict everything to the regions the second method was run on (for example its input peaks)
# Reports: fraction of bound sites and of unbound sites overlapped by the second set (the second must be clearly higher),
# enrichment ratio, and the reverse direction (fraction of the second set's calls that overlap a bound site).
# No universal pass cutoff exists; interpret the bound/unbound contrast, not a single percentage.
set -euo pipefail
[ $# -ge 3 ] || { echo "usage: $0 bound.bed unbound.bed other_calls.bed [regions.bed]" >&2; exit 2; }
BOUND=$1; UNBOUND=$2; OTHER=$3; REGIONS=${4:-}
for f in "$BOUND" "$UNBOUND" "$OTHER" ${REGIONS:+"$REGIONS"}; do [ -s "$f" ] || { echo "ERROR: missing or empty: $f" >&2; exit 2; }; done

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
prep() { cut -f1-3 "$1" | sort -u | { if [ -n "$REGIONS" ]; then bedtools intersect -u -a - -b "$REGIONS"; else cat; fi; }; }
prep "$BOUND" > "$T/bound.bed"; prep "$UNBOUND" > "$T/unbound.bed"; prep "$OTHER" > "$T/other.bed"

nb=$(wc -l < "$T/bound.bed"); nu=$(wc -l < "$T/unbound.bed"); no=$(wc -l < "$T/other.bed")
ob=$(bedtools intersect -u -a "$T/bound.bed" -b "$T/other.bed" | wc -l)
ou=$(bedtools intersect -u -a "$T/unbound.bed" -b "$T/other.bed" | wc -l)
oo=$(bedtools intersect -u -a "$T/other.bed" -b "$T/bound.bed" | wc -l)
awk -v nb="$nb" -v ob="$ob" -v nu="$nu" -v ou="$ou" -v no="$no" -v oo="$oo" 'BEGIN {
    fb = nb ? ob/nb : 0; fu = nu ? ou/nu : 0
    printf "bound sites overlapped by other set:   %d/%d = %.3f\n", ob, nb, fb
    printf "unbound sites overlapped by other set: %d/%d = %.3f\n", ou, nu, fu
    if (fu > 0) printf "enrichment (bound/unbound): %.1f\n", fb/fu; else print "enrichment (bound/unbound): undefined (unbound fraction 0)"
    printf "other-set calls overlapping a bound site: %d/%d = %.3f\n", oo, no, (no ? oo/no : 0) }'
