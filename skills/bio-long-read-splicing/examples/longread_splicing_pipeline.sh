#!/bin/bash
# Reference: checked 2026-09 on minimap2 2.31, samtools 1.24, bedtools 2.31.1, gffread 0.12.9, FLAIR 3.0.1, IsoQuant 4.0.0, SQANTI3 6.0.2 | Verify CLI flags if versions differ
# Long-read splicing analysis pipeline
#
# Workflow:
# 1. Splice-aware alignment with minimap2 (preset depends on platform)
# 2. Isoform discovery + quantification (IsoQuant, and FLAIR as the end-to-end alternative)
# 3. Classification + filtering with SQANTI3 (FSM/ISM/NIC/NNC + artifact flags)
#
# Set the platform with PLATFORM=hifi | ont | drna
#   hifi = PacBio HiFi/Iso-Seq;  ont = ONT cDNA (unstranded);  drna = ONT direct RNA (stranded)
# Inputs default to the names below; override any of them from the environment, e.g.
#   PLATFORM=ont REFERENCE=genome.fa GTF=anno.gtf FASTQ=reads.fq.gz SAMPLE=s1 bash longread_splicing_pipeline.sh
# Optional: SR_JUNCTIONS=SJ.out.tab (STAR short-read junctions) makes `flair correct` recover novel splice sites.

set -euo pipefail

PLATFORM=${PLATFORM:-hifi}
THREADS=${THREADS:-16}
REFERENCE=${REFERENCE:-reference.fa}
GTF=${GTF:-gencode.v45.annotation.gtf}
SAMPLE=${SAMPLE:-sample}
FASTQ=${FASTQ:-${SAMPLE}.fastq.gz}
SR_JUNCTIONS=${SR_JUNCTIONS:-}
OUTPUT_DIR=${OUTPUT_DIR:-longread_output_${SAMPLE}}

mkdir -p "${OUTPUT_DIR}"

# -uf only for reads that are all in transcript orientation (direct RNA, or orientation-fixed FLNC);
# on unstranded cDNA it makes ~half the reads align with false junctions.
if [ "${PLATFORM}" = "hifi" ]; then
    PRESET="splice:hq"
    DATA_TYPE="pacbio_ccs"
    STRANDED="none"
    JUNC_BONUS=9   # minimap2 default; higher values recode skipping reads of some microexons to inclusion on HiFi
elif [ "${PLATFORM}" = "ont" ]; then
    PRESET="splice -k14"
    DATA_TYPE="nanopore"
    STRANDED="none"
    JUNC_BONUS=16  # rescues annotated microexons below ~10 nt; see SKILL.md Microexons for why not higher
elif [ "${PLATFORM}" = "drna" ]; then
    PRESET="splice -uf -k14"
    DATA_TYPE="nanopore"
    STRANDED="forward"
    JUNC_BONUS=16
else
    echo "PLATFORM must be hifi, ont or drna" && exit 1
fi

# 1. Splice-aware alignment; --junc-bed makes minimap2 prefer annotated junctions (rescues annotated microexons; see SKILL.md Microexons for the bonus)
gffread "${GTF}" --bed -o "${OUTPUT_DIR}/annotation.bed12"

# $PRESET is intentionally unquoted: it expands to several minimap2 words
minimap2 -ax ${PRESET} \
    -t "${THREADS}" \
    --secondary=no \
    --junc-bed "${OUTPUT_DIR}/annotation.bed12" \
    --junc-bonus "${JUNC_BONUS}" \
    "${REFERENCE}" "${FASTQ}" | \
    samtools sort -@ "${THREADS}" -o "${OUTPUT_DIR}/${SAMPLE}_aligned.bam"
samtools index "${OUTPUT_DIR}/${SAMPLE}_aligned.bam"

# Orientation check (needs an alignment made WITHOUT -uf): ~1.0 = reads are in transcript orientation
# (-uf is then safe), ~0.5 = unoriented reads (never use -uf)
if [ "${PLATFORM}" != "drna" ]; then
    samtools view -F 2308 "${OUTPUT_DIR}/${SAMPLE}_aligned.bam" | \
        awk '{for(i=12;i<=NF;i++) if($i ~ /^ts:A:/){n++; if($i=="ts:A:+") p++}} END{if(n) printf "Spliced reads with ts:A:+ (transcript orientation): %.3f (n=%d)\n", p/n, n}'
fi

# Reads with an intron next to a soft clip: a --junc-bonus that is too high makes these (expect 0; they crash `flair correct`)
echo "Reads with an intron next to a soft clip: $(samtools view -F 2308 "${OUTPUT_DIR}/${SAMPLE}_aligned.bam" | awk '$6 ~ /N[0-9]+S$/ || $6 ~ /^[0-9]+S[0-9]+N/' | wc -l)"

# 2a. Isoform discovery and quantification with IsoQuant (entry point is `isoquant`; isoquant.py exists only in a git checkout)
isoquant \
    --reference "${REFERENCE}" \
    --genedb "${GTF}" \
    --bam "${OUTPUT_DIR}/${SAMPLE}_aligned.bam" \
    --data_type ${DATA_TYPE} \
    --stranded ${STRANDED} \
    --output "${OUTPUT_DIR}/isoquant" \
    --threads "${THREADS}" \
    --prefix "${SAMPLE}"
# --model_construction_strategy is auto-selected from --data_type; override only if needed
ISOQUANT_GTF=${OUTPUT_DIR}/isoquant/${SAMPLE}/${SAMPLE}.transcript_models.gtf
# SQANTI3 filter accepts transcript/exon rows only; IsoQuant also writes gene, CDS, UTR and start/stop_codon rows
awk -F'\t' '$3 == "transcript" || $3 == "exon"' "${ISOQUANT_GTF}" > "${OUTPUT_DIR}/isoquant_models.gtf"

# 2b. Alternative: FLAIR end-to-end pipeline
# Convert BAM to BED12 for FLAIR
bedtools bamtobed -bed12 -i "${OUTPUT_DIR}/${SAMPLE}_aligned.bam" > "${OUTPUT_DIR}/${SAMPLE}.bed"

# flair correct takes an annotation and/or orthogonal junctions (no --genome); novel splice sites need the junctions
SR_ARGS=()
if [ -n "${SR_JUNCTIONS}" ]; then SR_ARGS=(--junction_tab "${SR_JUNCTIONS}"); fi
flair correct \
    -q "${OUTPUT_DIR}/${SAMPLE}.bed" \
    -f "${GTF}" \
    "${SR_ARGS[@]}" \
    -o "${OUTPUT_DIR}/flair_corrected_${SAMPLE}" \
    -t "${THREADS}"

flair collapse \
    -q "${OUTPUT_DIR}/flair_corrected_${SAMPLE}_all_corrected.bed" \
    -r "${FASTQ}" \
    -g "${REFERENCE}" \
    -f "${GTF}" \
    -o "${OUTPUT_DIR}/flair_collapsed_${SAMPLE}" \
    -t "${THREADS}"

# 3. SQANTI3 classification (filter intra-priming and RT-switching flags before reporting)
# -o is a file prefix, -d the output directory; --report skip avoids the slow PDF/HTML report
sqanti3_qc.py \
    --isoforms "${OUTPUT_DIR}/isoquant_models.gtf" \
    --refGTF "${GTF}" \
    --refFasta "${REFERENCE}" \
    -o sqanti3 -d "${OUTPUT_DIR}/sqanti3" \
    --aligner_choice minimap2 \
    --report skip \
    --cpus "${THREADS}"

# --skip_report: the filter's R report step fails on small inputs
sqanti3_filter.py rules \
    --sqanti_class "${OUTPUT_DIR}/sqanti3/sqanti3_classification.txt" \
    --filter_gtf "${OUTPUT_DIR}/isoquant_models.gtf" \
    -o sqanti3_filtered -d "${OUTPUT_DIR}/sqanti3_filtered" \
    --skip_report

echo "Pipeline complete. Outputs in ${OUTPUT_DIR}/"
echo "  IsoQuant transcript models: ${ISOQUANT_GTF}"
echo "  FLAIR isoforms: ${OUTPUT_DIR}/flair_collapsed_${SAMPLE}.isoforms.fa"
echo "  SQANTI3 classification: ${OUTPUT_DIR}/sqanti3/sqanti3_classification.txt"
echo "  SQANTI3 filtered (rt-switching, intra-priming removed): ${OUTPUT_DIR}/sqanti3_filtered/sqanti3_filtered.filtered.gtf"
