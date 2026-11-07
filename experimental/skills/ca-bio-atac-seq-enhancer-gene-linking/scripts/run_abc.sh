#!/bin/bash
# Standalone ABC run: accessibility (+ H3K27ac, + Hi-C) -> unthresholded and thresholded enhancer-gene links.
# Mirrors the ABC Snakemake rules (call_macs_peaks, make_candidate_regions, create_neighborhoods,
# create_predictions, filter_predictions) with default parameters. Tested against ABC v1.1.2 and main (92ac503).
# Tools on PATH: macs2 (only when PEAKS is unset), bedtools, python with the ABC environment.
#
# Required env vars:
#   ABC_REPO      clone of broadinstitute/ABC-Enhancer-Gene-Prediction (scripts in workflow/scripts/)
#   ACCESS_BAM    accessibility reads: DNase BAM, ATAC BAM, or tagAlign for ATAC (comma-separate replicates)
#   OUTDIR        output directory
# Optional env vars (defaults in brackets):
#   ACCESS_TYPE   DHS | ATAC [ATAC]; sets --DHS/--ATAC, --accessibility_feature and the threshold row
#   H3K27AC_BAM   H3K27ac ChIP BAM; unset runs the accessibility-only model
#   HIC_TYPE      none | hic | avg | juicebox | bedpe [none]; none = powerlaw contact
#   HIC_FILE      .hic file or URL (hic), or directory (avg/juicebox/bedpe); required unless HIC_TYPE=none
#   HIC_RES       Hi-C resolution in bp [5000]
#   PEAKS         MACS2 narrowPeak made with --call-summits; unset runs macs2 as the ABC pipeline does
#   CELL_TYPE     label written to outputs [sample]
#   CHROM_SIZES   chrom sizes TSV [ABC_REPO/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv]
#   GENES_BED     BED6 gene bounds [ABC_REPO/reference/hg38/CollapsedGeneBounds.hg38.bed]
#   TSS_BED       500 bp TSS regions, forced into candidates [ABC_REPO/reference/hg38/CollapsedGeneBounds.hg38.TSS500bp.bed]
#   BLOCKLIST     regions removed from candidates [ABC_REPO/reference/hg38/GRCh38_unified_blacklist.bed]
#   UBIQ_GENES    ubiquitously expressed genes [ABC_REPO/reference/UbiquitouslyExpressedGenes.txt]
#   QNORM         quantile-normalisation reference, or "none" [ABC_REPO/reference/EnhancersQNormRef.K562.txt]
#   MACS_GENOME   macs2 -g value [hs]
#   THRESHOLD     score cut-off; unset looks it up in ABC_REPO/reference/abc_thresholds.tsv by ACCESS_TYPE, H3K27ac and Hi-C type
#   HIC_GAMMA, HIC_SCALE  powerlaw parameters [ABC config average-Hi-C values 1.024238616787792, 5.9594510043736655]

set -euo pipefail

: "${ABC_REPO:?set ABC_REPO}" "${ACCESS_BAM:?set ACCESS_BAM}" "${OUTDIR:?set OUTDIR}"
ACCESS_TYPE=${ACCESS_TYPE:-ATAC}
H3K27AC_BAM=${H3K27AC_BAM:-}
HIC_TYPE=${HIC_TYPE:-none}
HIC_FILE=${HIC_FILE:-}
HIC_RES=${HIC_RES:-5000}
PEAKS=${PEAKS:-}
CELL_TYPE=${CELL_TYPE:-sample}
REF=$ABC_REPO/reference
CHROM_SIZES=${CHROM_SIZES:-$REF/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv}
GENES_BED=${GENES_BED:-$REF/hg38/CollapsedGeneBounds.hg38.bed}
TSS_BED=${TSS_BED:-$REF/hg38/CollapsedGeneBounds.hg38.TSS500bp.bed}
BLOCKLIST=${BLOCKLIST:-$REF/hg38/GRCh38_unified_blacklist.bed}
UBIQ_GENES=${UBIQ_GENES:-$REF/UbiquitouslyExpressedGenes.txt}
QNORM=${QNORM:-$REF/EnhancersQNormRef.K562.txt}
MACS_GENOME=${MACS_GENOME:-hs}
HIC_GAMMA=${HIC_GAMMA:-1.024238616787792}
HIC_SCALE=${HIC_SCALE:-5.9594510043736655}
SCRIPTS=$ABC_REPO/workflow/scripts

[[ -f $SCRIPTS/predict.py ]] || { echo "no $SCRIPTS/predict.py: ABC_REPO must be a v1.0+ or main clone" >&2; exit 2; }
case $ACCESS_TYPE in DHS|ATAC) ;; *) echo "ACCESS_TYPE must be DHS or ATAC" >&2; exit 2;; esac
case $HIC_TYPE in
    none) [[ -z $HIC_FILE ]] || { echo "HIC_FILE is set but HIC_TYPE=none" >&2; exit 2; } ;;
    hic|avg|juicebox|bedpe) : "${HIC_FILE:?HIC_TYPE=$HIC_TYPE needs HIC_FILE}" ;;
    *) echo "HIC_TYPE must be none, hic, avg, juicebox or bedpe" >&2; exit 2 ;;
esac

mkdir -p "$OUTDIR"/{peaks,neighborhoods,predictions}
ACCESS_FIRST=${ACCESS_BAM%%,*}
IFS=, read -ra ACCESS_LIST <<< "$ACCESS_BAM"

# 1. Chromosome BED, then peaks (ABC's MACS2 settings; --call-summits is required by candidate-region step)
awk 'BEGIN{OFS="\t"} NF>0{print $1,0,$2}' "$CHROM_SIZES" > "$OUTDIR/peaks/chrom_sizes.bed"
if [[ -z $PEAKS ]]; then
    FMT=AUTO; [[ $ACCESS_FIRST == *tagAlign* ]] && FMT=BED
    macs2 callpeak -f "$FMT" -g "$MACS_GENOME" -p 0.1 -n macs2 --shift -75 --extsize 150 --nomodel \
        --keep-dup all --call-summits --outdir "$OUTDIR/peaks" -t "${ACCESS_LIST[@]}"
    PEAKS=$OUTDIR/peaks/macs2_peaks.narrowPeak
fi
bedtools intersect -u -a "$PEAKS" -b "$OUTDIR/peaks/chrom_sizes.bed" \
    | bedtools sort -faidx "$CHROM_SIZES" -i stdin > "$OUTDIR/peaks/peaks.sorted.narrowPeak"

# 2. Candidate regions: summit +/- 250 bp, top 150000 peaks, blocklist removed, TSS regions added
python "$SCRIPTS/makeCandidateRegions.py" \
    --narrowPeak "$OUTDIR/peaks/peaks.sorted.narrowPeak" \
    --accessibility "${ACCESS_LIST[@]}" \
    --outDir "$OUTDIR/peaks" \
    --chrom_sizes "$CHROM_SIZES" --chrom_sizes_bed "$OUTDIR/peaks/chrom_sizes.bed" \
    --regions_blocklist "$BLOCKLIST" --regions_includelist "$TSS_BED" \
    --peakExtendFromSummit 250 --nStrongestPeak 150000
CANDIDATES=$OUTDIR/peaks/peaks.sorted.narrowPeak.candidateRegions.bed

# 3. Neighborhoods: Activity from BAM read counts, quantile-normalised
bedtools intersect -u -a "$GENES_BED" -b "$OUTDIR/peaks/chrom_sizes.bed" \
    | bedtools sort -faidx "$CHROM_SIZES" -i stdin | uniq > "$OUTDIR/neighborhoods/genes.processed.bed"
NB_ARGS=()
[[ $QNORM == none ]] || NB_ARGS+=(--qnorm "$QNORM")
[[ -z $H3K27AC_BAM ]] || NB_ARGS+=(--H3K27ac "$H3K27AC_BAM")
python "$SCRIPTS/run.neighborhoods.py" \
    --candidate_enhancer_regions "$CANDIDATES" \
    "--$ACCESS_TYPE" "$ACCESS_BAM" --default_accessibility_feature "$ACCESS_TYPE" \
    --genes "$OUTDIR/neighborhoods/genes.processed.bed" \
    --ubiquitously_expressed_genes "$UBIQ_GENES" \
    --chrom_sizes "$CHROM_SIZES" --chrom_sizes_bed "$OUTDIR/peaks/chrom_sizes.bed" \
    --cellType "$CELL_TYPE" --outdir "$OUTDIR/neighborhoods" "${NB_ARGS[@]}"

# 4. Predictions: Activity x Contact for all pairs within 5 Mb
SCORE=ABC.Score; HIC_ARGS=(); THR_TYPE=powerlaw
if [[ $HIC_TYPE != none ]]; then
    HIC_ARGS=(--hic_file "$HIC_FILE" --hic_type "$HIC_TYPE" --hic_resolution "$HIC_RES" --scale_hic_using_powerlaw)
    THR_TYPE=intact_hic; [[ $HIC_TYPE == avg ]] && THR_TYPE=avg
else
    SCORE=powerlaw.Score
fi
python "$SCRIPTS/predict.py" \
    --enhancers "$OUTDIR/neighborhoods/EnhancerList.txt" --genes "$OUTDIR/neighborhoods/GeneList.txt" \
    --chrom_sizes "$CHROM_SIZES" --accessibility_feature "$ACCESS_TYPE" --cellType "$CELL_TYPE" \
    --score_column "$SCORE" --hic_gamma "$HIC_GAMMA" --hic_scale "$HIC_SCALE" \
    --hic_pseudocount_distance 5000 "${HIC_ARGS[@]}" --outdir "$OUTDIR/predictions"

# 5. Threshold: calibrated per input combination, then drop non-self promoters
if [[ -z ${THRESHOLD:-} ]]; then
    HAS_H3=FALSE; [[ -n $H3K27AC_BAM ]] && HAS_H3=TRUE
    THRESHOLD=$(awk -F'\t' -v a="$ACCESS_TYPE" -v h="$HAS_H3" -v t="$THR_TYPE" \
        '$1==a && $2==h && $3==t {print $4}' "$REF/abc_thresholds.tsv")
    [[ -n $THRESHOLD ]] || THRESHOLD=0.02
fi
echo "score column $SCORE, threshold $THRESHOLD"
python "$SCRIPTS/filter_predictions.py" \
    --output_tsv_file "$OUTDIR/predictions/EnhancerPredictionsFull_threshold$THRESHOLD.tsv" \
    --output_slim_tsv_file "$OUTDIR/predictions/EnhancerPredictions_threshold$THRESHOLD.tsv" \
    --output_bed_file "$OUTDIR/predictions/EnhancerPredictionsFull_threshold$THRESHOLD.bedpe.gz" \
    --output_gene_stats_file "$OUTDIR/predictions/GenePredictionStats_threshold$THRESHOLD.tsv" \
    --pred_file "$OUTDIR/predictions/EnhancerPredictionsAllPutative.tsv.gz" \
    --pred_nonexpressed_file "$OUTDIR/predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz" \
    --score_column "$SCORE" --threshold "$THRESHOLD" \
    --include_self_promoter True --only_expressed_genes False
