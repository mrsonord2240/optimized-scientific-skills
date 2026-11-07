#!/usr/bin/env bash
# Tested contract: WASP d3b8447, GATK 4.6, samtools/bcftools 1.21, bowtie2 2.5.
# WASP correction + GATK ASEReadCounter + phase-oriented peak aggregation.

set -euo pipefail

if [[ $# -ne 7 ]]; then
    echo "usage: $0 BAM VCF GENOME_FASTA BOWTIE2_INDEX PEAKS_BED WASP_DIR OUTPUT_DIR" >&2
    exit 64
fi

BAM=$1
VCF=$2
GENOME_FA=$3
BT2_INDEX=$4
PEAKS=$5
WASP_DIR=$6
OUTDIR=$7
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

for path in "$BAM" "$VCF" "$GENOME_FA" "${GENOME_FA}.fai" "$PEAKS"; do
    if [[ ! -r "$path" ]]; then
        echo "required input is not readable: $path" >&2
        exit 66
    fi
done
for executable in samtools bcftools bgzip tabix snp2h5 bowtie2 gatk python python3; do
    if ! command -v "$executable" >/dev/null 2>&1; then
        echo "required executable is unavailable: $executable" >&2
        exit 69
    fi
done
if [[ ! -r "$WASP_DIR/mapping/find_intersecting_snps.py" || ! -r "$WASP_DIR/mapping/filter_remapped_reads.py" ]]; then
    echo "WASP_DIR does not contain the required mapping scripts: $WASP_DIR" >&2
    exit 66
fi
if [[ -e "$OUTDIR" ]]; then
    echo "output directory already exists; choose a new path or archive/remove it explicitly: $OUTDIR" >&2
    exit 73
fi

samtools quickcheck "$BAM"
BAM_HEADER=$(samtools view -H "$BAM")
mapfile -t BAM_SAMPLES < <(
    awk -F '\t' '
        $1 == "@RG" {
            for (i = 2; i <= NF; i++) {
                if ($i ~ /^SM:/) {
                    sub(/^SM:/, "", $i)
                    if ($i != "") print $i
                }
            }
        }
    ' <<< "$BAM_HEADER" | LC_ALL=C sort -u
)
if [[ ${#BAM_SAMPLES[@]} -ne 1 ]]; then
    echo "BAM must contain read groups for exactly one non-empty SM sample; found ${#BAM_SAMPLES[@]}" >&2
    exit 65
fi
SAMPLE=${BAM_SAMPLES[0]}
if [[ ! "$SAMPLE" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "BAM SM contains characters unsafe for output filenames: $SAMPLE" >&2
    exit 65
fi
mapfile -t VCF_SAMPLES < <(bcftools query -l "$VCF")
MATCHES=0
for vcf_sample in "${VCF_SAMPLES[@]}"; do
    [[ "$vcf_sample" == "$SAMPLE" ]] && MATCHES=$((MATCHES + 1))
done
if [[ $MATCHES -ne 1 ]]; then
    echo "BAM sample $SAMPLE must match exactly one VCF sample; matches=$MATCHES" >&2
    exit 65
fi

OUT_PARENT=$(dirname "$OUTDIR")
mkdir -p "$OUT_PARENT"
STAGE=$(mktemp -d "${OUTDIR}.tmp.XXXXXX")
COMPLETED=0
cleanup() {
    if [[ $COMPLETED -eq 0 && -d "$STAGE" ]]; then
        rm -rf -- "$STAGE"
    fi
}
trap cleanup EXIT
mkdir -p "$STAGE/wasp/snp_h5" "$STAGE/ase" "$STAGE/peak_aggregation"

# Select the BAM sample before applying the heterozygous filter. This prevents
# another sample's heterozygous genotype from admitting a target-homozygous site.
TARGET_VCF="$STAGE/wasp/${SAMPLE}.target.vcf.gz"
FILTERED_VCF="$STAGE/wasp/${SAMPLE}.phased_het.vcf.gz"
bcftools view --samples "$SAMPLE" -Oz -o "$TARGET_VCF" "$VCF"
tabix -f -p vcf "$TARGET_VCF"
bcftools view -m2 -M2 -v snps -i 'GT="het"' \
    -Oz -o "$FILTERED_VCF" "$TARGET_VCF"
tabix -f -p vcf "$FILTERED_VCF"
if [[ $(bcftools query -l "$FILTERED_VCF") != "$SAMPLE" ]]; then
    echo "filtered VCF did not retain exactly the target sample $SAMPLE" >&2
    exit 65
fi
if [[ $(bcftools index -n "$FILTERED_VCF") -eq 0 ]]; then
    echo "target sample has no biallelic heterozygous SNPs after filtering" >&2
    exit 65
fi

# Build an explicit, auditable phase-orientation table. Missing PS or unphased
# GT is a hard stop: those sites cannot be pooled into a shared haplotype.
HAPLOTYPE_MAP="$STAGE/ase/${SAMPLE}.haplotype_map.tsv"
printf 'contig\tposition\trefAllele\taltAllele\tphaseSet\thaplotype1Allele\n' > "$HAPLOTYPE_MAP"
bcftools query --samples "$SAMPLE" \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%PS]\n' "$FILTERED_VCF" | \
    awk -F '\t' -v OFS='\t' '
        $5 != "0|1" && $5 != "1|0" {
            print "unphased or unsupported target genotype at " $1 ":" $2 ": " $5 > "/dev/stderr"
            bad = 1
            next
        }
        $6 == "" || $6 == "." {
            print "missing phase-set PS at " $1 ":" $2 > "/dev/stderr"
            bad = 1
            next
        }
        { print $1, $2, $3, $4, $1 ":" $6, ($5 == "0|1" ? "REF" : "ALT") }
        END { if (bad) exit 65 }
    ' >> "$HAPLOTYPE_MAP"
if [[ $(wc -l < "$HAPLOTYPE_MAP") -le 1 ]]; then
    echo "haplotype map contains no phased target variants" >&2
    exit 65
fi

cut -f1,2 "${GENOME_FA}.fai" > "$STAGE/wasp/chromInfo.txt"
snp2h5 --chrom "$STAGE/wasp/chromInfo.txt" \
    --format vcf \
    --snp_tab "$STAGE/wasp/snp_h5/snp_tab.h5" \
    --snp_index "$STAGE/wasp/snp_h5/snp_index.h5" \
    --haplotype "$STAGE/wasp/snp_h5/haps.h5" \
    "$FILTERED_VCF"

python "$WASP_DIR/mapping/find_intersecting_snps.py" \
    --is_paired_end \
    --is_sorted \
    --output_dir "$STAGE/wasp" \
    --snp_tab "$STAGE/wasp/snp_h5/snp_tab.h5" \
    --snp_index "$STAGE/wasp/snp_h5/snp_index.h5" \
    --haplotype "$STAGE/wasp/snp_h5/haps.h5" \
    --samples "$SAMPLE" \
    "$BAM"

bowtie2 -x "$BT2_INDEX" \
    -1 "$STAGE/wasp/${SAMPLE}.remap.fq1.gz" \
    -2 "$STAGE/wasp/${SAMPLE}.remap.fq2.gz" \
    -p "${ASE_THREADS:-8}" -S "$STAGE/wasp/${SAMPLE}.remap.sam"
samtools view -bS "$STAGE/wasp/${SAMPLE}.remap.sam" | \
    samtools sort -o "$STAGE/wasp/${SAMPLE}.remap.bam"
samtools index "$STAGE/wasp/${SAMPLE}.remap.bam"

python "$WASP_DIR/mapping/filter_remapped_reads.py" \
    "$STAGE/wasp/${SAMPLE}.to.remap.bam" \
    "$STAGE/wasp/${SAMPLE}.remap.bam" \
    "$STAGE/wasp/${SAMPLE}.kept.bam"

# A merge of individually sorted inputs is not guaranteed to remain sorted.
# Sort after the merge, add a read group only to orphan reads, then validate the
# final BAM's order, sample identity, and index.
samtools merge -f "$STAGE/wasp/${SAMPLE}.merged.unsorted.bam" \
    "$STAGE/wasp/${SAMPLE}.kept.bam" \
    "$STAGE/wasp/${SAMPLE}.keep.bam"
samtools sort -o "$STAGE/wasp/${SAMPLE}.merged.sorted.bam" \
    "$STAGE/wasp/${SAMPLE}.merged.unsorted.bam"
samtools addreplacerg \
    -r "ID:WASP_${SAMPLE}" -r "SM:${SAMPLE}" -r 'PL:ILLUMINA' \
    -m orphan_only \
    -o "$STAGE/wasp/${SAMPLE}.wasp.bam" \
    "$STAGE/wasp/${SAMPLE}.merged.sorted.bam"
samtools quickcheck "$STAGE/wasp/${SAMPLE}.wasp.bam"
FINAL_SORT_ORDER=$(samtools view -H "$STAGE/wasp/${SAMPLE}.wasp.bam" | awk -F '\t' '
    $1 == "@HD" { for (i = 2; i <= NF; i++) if ($i ~ /^SO:/) { sub(/^SO:/, "", $i); print $i } }
')
if [[ "$FINAL_SORT_ORDER" != "coordinate" ]]; then
    echo "final WASP BAM is not coordinate sorted: SO=$FINAL_SORT_ORDER" >&2
    exit 70
fi
mapfile -t FINAL_SAMPLES < <(samtools view -H "$STAGE/wasp/${SAMPLE}.wasp.bam" | awk -F '\t' '
    $1 == "@RG" { for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) { sub(/^SM:/, "", $i); print $i } }
' | LC_ALL=C sort -u)
if [[ ${#FINAL_SAMPLES[@]} -ne 1 || "${FINAL_SAMPLES[0]}" != "$SAMPLE" ]]; then
    echo "final WASP BAM read groups do not resolve uniquely to $SAMPLE" >&2
    exit 70
fi
samtools index "$STAGE/wasp/${SAMPLE}.wasp.bam"

ASE_COUNTS="$STAGE/ase/${SAMPLE}.ase_counts.tsv"
gatk ASEReadCounter \
    -I "$STAGE/wasp/${SAMPLE}.wasp.bam" \
    -V "$FILTERED_VCF" \
    -R "$GENOME_FA" \
    -O "$ASE_COUNTS" \
    --output-format TABLE \
    --min-mapping-quality 30 \
    --min-base-quality 20
if ! awk 'NR > 1 && NF > 0 { found = 1 } END { exit(found ? 0 : 1) }' "$ASE_COUNTS"; then
    echo "GATK ASEReadCounter wrote no data rows; refusing to aggregate" >&2
    exit 70
fi

python3 "$SCRIPT_DIR/aggregate_peak_ase.py" \
    --ase-counts "$ASE_COUNTS" \
    --peaks "$PEAKS" \
    --haplotype-map "$HAPLOTYPE_MAP" \
    --output "$STAGE/peak_aggregation/${SAMPLE}.peak_ase.tsv"

mv -- "$STAGE" "$OUTDIR"
COMPLETED=1
trap - EXIT
echo "Done. Validated outputs: $OUTDIR/{wasp,ase,peak_aggregation}/"
