#!/usr/bin/env bash
# Build the leafviz annotation and results .RData from leafcutter output, then start the Shiny app.
# Usage: leafviz_run.sh LEAFCUTTER_DIR ANNOTATION.gtf GROUPS.txt PERIND_NUMERS.counts.gz DS_PREFIX OUT.RData
#   LEAFCUTTER_DIR  the leafcutter repo clone (contains leafviz/)   GROUPS.txt  the support file given to leafcutter_ds.R (sample<TAB>condition)
#   DS_PREFIX       the -o prefix of leafcutter_ds.R (reads DS_PREFIX_cluster_significance.txt and DS_PREFIX_effect_sizes.txt)
#   Env: LAUNCH=0 builds the files and stops before the app. Needs perl and Rscript on PATH.
# annotation_code = prefix "annot" (_all_exons.txt.gz, _all_introns.bed.gz, _fiveprime.bed.gz, _threeprime.bed.gz) written next to OUT.RData;
# build it from the GTF version used in the differential analysis.
set -euo pipefail
[ $# -eq 6 ] || { sed -n '2,7p' "$0" >&2; exit 2; }
lc=$(cd "$1" && pwd); gtf=$2; groups=$3; counts=$4; ds=$5; rdata=$6
out=$(cd "$(dirname "$rdata")" && pwd)/$(basename "$rdata")
annot=$(dirname "$out")/annot

perl "$lc/leafviz/gtf2leafcutter.pl" -o "$annot" "$gtf"
Rscript "$lc/leafviz/prepare_results.R" -o "$out" -m "$groups" \
    "$counts" "${ds}_cluster_significance.txt" "${ds}_effect_sizes.txt" "$annot"
[ -s "$out" ] || { echo "prepare_results.R wrote no $out" >&2; exit 1; }
[ "${LAUNCH:-1}" = 1 ] || { echo "OK: $out built; app not started (LAUNCH=0)"; exit 0; }

# runApp() uses the working directory: start from leafviz/, pass the .RData by absolute path
cd "$lc/leafviz" && exec Rscript run_leafviz.R "$out"    # prints "Listening on http://127.0.0.1:<port>"
