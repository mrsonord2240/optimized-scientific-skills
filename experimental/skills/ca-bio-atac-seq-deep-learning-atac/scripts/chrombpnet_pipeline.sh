#!/bin/bash
# chromBPNet 1.0.1 (env `chrombpnet`, see SKILL.md): splits -> nonpeaks -> bias model -> accessibility model -> QC -> optional variant scoring.
# CPU-only (TF 2.8). Uses `bias train` / `train`, which stop after training, so the slow DeepSHAP interpret tail of
# `bias pipeline` / `pipeline` is skipped; run interpretation separately (SKILL.md step 6).
#
# usage: chrombpnet_pipeline.sh [BAM] [PEAKS] [GENOME] [SIZES] [OUTDIR] [ATAC|DNASE] [BIAS_THRESH]
#   PEAKS      10-column narrowPeak (summit in column 10)
#   Env: TEST_CHROMS ("chr1 chr3 chr6"), VALID_CHROMS ("chr8 chr20"); everything else trains.
#        TRAIN_ARGS  extra flags for both train steps, e.g. "-e 2 -es 1" for a smoke run
#        VARIANTS    5-column headerless list (chr pos ref alt variant_id); needs VARIANT_SCORER=<variant-scorer clone>
#        ALLOW_FAILED_QC=1 continues to variant scoring although the QC gate failed
#        RESUME=1    skip every step whose output already exists (delete a half-written bias/ or model/ directory first)

set -euo pipefail

BAM=${1:-atac.dedup.bam}
PEAKS=${2:-peaks.narrowPeak}
GENOME=${3:-hg38.fa}
SIZES=${4:-hg38.chrom.sizes}
OUTDIR=${5:-chrombpnet_out}
DATA_TYPE=${6:-ATAC}
BIAS_THRESH=${7:-0.5}
TEST_CHROMS=${TEST_CHROMS:-chr1 chr3 chr6}
VALID_CHROMS=${VALID_CHROMS:-chr8 chr20}
TRAIN_ARGS=${TRAIN_ARGS:-}
RESUME=${RESUME:-0}
skip() { [ "$RESUME" = 1 ] && [ -s "$1" ]; }
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# 0. Fail fast on missing inputs, before any long step
for t in chrombpnet bedtools bedGraphToBigWig; do
    command -v "$t" >/dev/null || { echo "$t not on PATH (activate the chrombpnet env)" >&2; exit 1; }
done
for f in "$BAM" "$PEAKS" "$GENOME" "$SIZES"; do
    [ -s "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done
[ "$(head -n1 "$PEAKS" | awk -F'\t' '{print NF}')" -ge 10 ] || { echo "$PEAKS must be 10-column narrowPeak" >&2; exit 1; }
if [ -n "${VARIANTS:-}" ]; then
    [ -s "$VARIANTS" ] || { echo "missing VARIANTS file: $VARIANTS" >&2; exit 1; }
    [ "$(head -n1 "$VARIANTS" | awk -F'\t' '{print NF}')" -eq 5 ] || { echo "$VARIANTS needs 5 tab-separated columns: chr pos ref alt variant_id" >&2; exit 1; }
    [ -f "${VARIANT_SCORER:-}/src/variant_scoring.py" ] || { echo "set VARIANT_SCORER to a kundajelab/variant-scorer clone" >&2; exit 1; }
fi
[ "$RESUME" = 1 ] || [ ! -e "$OUTDIR/nonpeaks_auxiliary" ] || { echo "$OUTDIR/nonpeaks_auxiliary exists; prep nonpeaks refuses to overwrite" >&2; exit 1; }
[ "$RESUME" = 1 ] || { [ ! -e "$OUTDIR/bias/models" ] && [ ! -e "$OUTDIR/model/models" ]; } || { echo "$OUTDIR/{bias,model} already trained; use a fresh OUTDIR" >&2; exit 1; }

mkdir -p "$OUTDIR"/{bias,model,splits,preds,variants}

# 1. Chromosome splits (JSON); -tcr = test, -vcr = validation, all remaining chromosomes train
# shellcheck disable=SC2086
skip "$OUTDIR/splits/fold_0.json" || chrombpnet prep splits -c "$SIZES" -tcr $TEST_CHROMS -vcr $VALID_CHROMS -op "$OUTDIR/splits/fold_0"

# 2. GC-matched non-peak background -> $OUTDIR/nonpeaks_negatives.bed
skip "$OUTDIR/nonpeaks_negatives.bed" || chrombpnet prep nonpeaks -g "$GENOME" -c "$SIZES" -p "$PEAKS" -fl "$OUTDIR/splits/fold_0.json" -o "$OUTDIR/nonpeaks"
NONPEAKS=$OUTDIR/nonpeaks_negatives.bed

# 3. Bias model from non-peak regions (same library; not a separate naked-DNA control)
# shellcheck disable=SC2086
skip "$OUTDIR/bias/models/bias.h5" || chrombpnet bias train -ibam "$BAM" -d "$DATA_TYPE" -g "$GENOME" -c "$SIZES" -p "$PEAKS" -n "$NONPEAKS" \
    -fl "$OUTDIR/splits/fold_0.json" -b "$BIAS_THRESH" $TRAIN_ARGS -o "$OUTDIR/bias/"

# 4. Accessibility model with bias correction -> model/models/{chrombpnet,chrombpnet_nobias}.h5
# shellcheck disable=SC2086
skip "$OUTDIR/model/models/chrombpnet_nobias.h5" || chrombpnet train -ibam "$BAM" -d "$DATA_TYPE" -g "$GENOME" -c "$SIZES" -p "$PEAKS" -n "$NONPEAKS" \
    -fl "$OUTDIR/splits/fold_0.json" -b "$OUTDIR/bias/models/bias.h5" $TRAIN_ARGS -o "$OUTDIR/model/"

# 5. QC on held-out test-chromosome peaks: corrected-model metrics (-bw observed) and Tn5 bias response
awk -F'\t' -v OFS='\t' -v want="$TEST_CHROMS" 'BEGIN{n=split(want,a," ");for(i=1;i<=n;i++)t[a[i]]=1} ($1 in t)' "$PEAKS" > "$OUTDIR/preds/test_peaks.narrowPeak"
[ -s "$OUTDIR/preds/test_peaks.narrowPeak" ] || { echo "no peaks on test chromosomes: $TEST_CHROMS" >&2; exit 1; }
skip "$OUTDIR/preds/qc_chrombpnet_metrics.json" || chrombpnet pred_bw -cm "$OUTDIR/model/models/chrombpnet.h5" -cmb "$OUTDIR/model/models/chrombpnet_nobias.h5" \
    -bm "$OUTDIR/bias/models/bias.h5" -r "$OUTDIR/preds/test_peaks.narrowPeak" -g "$GENOME" -c "$SIZES" \
    -bw "$OUTDIR/model/auxiliary/data_unstranded.bw" -op "$OUTDIR/preds/qc"
MOTIFS=$(DT="$DATA_TYPE" python -c 'import chrombpnet,os;print(os.path.join(os.path.dirname(chrombpnet.__file__),"data","motif_to_pwm.%s.tsv"%os.environ["DT"]))')
# footprints inserts 6 sequences into every test-chromosome non-peak window; cap at ~1000 windows to bound CPU time
awk -F'\t' -v want="$TEST_CHROMS" 'BEGIN{n=split(want,a," ");for(i=1;i<=n;i++)t[a[i]]=1} ($1 in t)' "$NONPEAKS" > "$OUTDIR/preds/test_nonpeaks.all.bed"
awk -v n="$(wc -l < "$OUTDIR/preds/test_nonpeaks.all.bed")" 'BEGIN{k=int(n/1000)+1} (NR-1)%k==0' "$OUTDIR/preds/test_nonpeaks.all.bed" > "$OUTDIR/preds/test_nonpeaks.bed"
skip "$OUTDIR/preds/nobias_max_bias_response.txt" || chrombpnet footprints -m "$OUTDIR/model/models/chrombpnet_nobias.h5" -r "$OUTDIR/preds/test_nonpeaks.bed" -g "$GENOME" \
    -fl "$OUTDIR/splits/fold_0.json" -op "$OUTDIR/preds/nobias" -pwm_f "$MOTIFS"

QC_OK=1
python "$HERE/chrombpnet_qc_gate.py" "$OUTDIR/model/evaluation/bias_metrics.json" \
    "$OUTDIR/preds/qc_chrombpnet_metrics.json" "$OUTDIR/preds/nobias_max_bias_response.txt" || QC_OK=0
[ "$QC_OK" = 1 ] || [ "${ALLOW_FAILED_QC:-0}" = 1 ] || { echo "QC gate failed: do not use this model for variants or motifs" >&2; exit 1; }

# 6. Optional variant scoring with the separate kundajelab/variant-scorer repo
if [ -n "${VARIANTS:-}" ]; then
    python "$VARIANT_SCORER/src/variant_scoring.py" \
        --model "$OUTDIR/model/models/chrombpnet_nobias.h5" \
        --list "$VARIANTS" --genome "$GENOME" --chrom_sizes "$SIZES" \
        --out_prefix "$OUTDIR/variants/predictions" --no_hdf5
    echo "Variant scores: $OUTDIR/variants/predictions.variant_scores.tsv (logfc = log2(alt/ref counts))"
fi

echo "Models: $OUTDIR/model/models/chrombpnet_nobias.h5  bias: $OUTDIR/bias/models/bias.h5"
