#!/bin/bash
# GATK-style hard filters on a germline VCF: split SNPs and indels, filter each with its own thresholds,
# keep sites whose RankSum annotations are missing (the GATK "missing => PASS" rule, so hom-alt sites survive),
# merge, and print Ti/Tv before and after.
#
# Usage:  [ASSAY=wgs|wes] [ALLOW_DL_CALLER=1] [ALLOW_DROP=1] hard_filter.sh <in.vcf[.gz]> [out_prefix=filtered]
# Output: <prefix>_all_filtered.vcf.gz (+ .tbi), plus <prefix>_snps_filtered.vcf.gz and <prefix>_indels_filtered.vcf.gz
# Needs:  bcftools. Checked on bcftools 1.24.
# Thresholds (GATK lenient starting points; a record fails if ANY applies):
#   SNP    QUAL<30  QD<2  FS>60   MQ<40  MQRankSum<-12.5  ReadPosRankSum<-8   SOR>3
#   indel  QUAL<30  QD<2  FS>200  ReadPosRankSum<-20  SOR>10          (MQ and MQRankSum are not used for indels)
# Traps (exit 2):
#   DeepVariant or DRAGEN header     their calibrated fields already encode quality; ALLOW_DL_CALLER=1 overrides
#   a GATK annotation absent from the INFO header (the thresholds mean nothing without it)
#   SNP + indel records != all records   MNPs, mixed or multiallelic SNP+indel records would be dropped or doubled;
#                                         normalise first (bcftools norm -m- --atomize) or ALLOW_DROP=1
#   empty output                     exit 3

set -euo pipefail

IN=${1:-}
PREFIX=${2:-filtered}
if [ -z "$IN" ]; then
    echo "Usage: [ASSAY=wgs|wes] hard_filter.sh <in.vcf[.gz]> [out_prefix]" >&2
    exit 1
fi

HDR=$(bcftools view -h "$IN")
if grep -Eiq 'deepvariant|dragen' <<<"$HDR" && [ "${ALLOW_DL_CALLER:-0}" != 1 ]; then
    echo "ERROR: $IN looks like DeepVariant or DRAGEN output. Filter on the caller's own QUAL and FILTER, not GATK hard filters (ALLOW_DL_CALLER=1 overrides)." >&2
    exit 2
fi
for tag in QD FS SOR MQ MQRankSum ReadPosRankSum; do
    if ! grep -q "^##INFO=<ID=$tag," <<<"$HDR"; then
        echo "ERROR: INFO/$tag is not in the header of $IN; GATK hard filters need the HaplotypeCaller annotations." >&2
        exit 2
    fi
done

bcftools view -v snps "$IN" -Oz -o "${PREFIX}_snps_raw.vcf.gz"
bcftools view -v indels "$IN" -Oz -o "${PREFIX}_indels_raw.vcf.gz"

n_total=$(bcftools view -H "$IN" | wc -l)
n_snp=$(bcftools view -H "${PREFIX}_snps_raw.vcf.gz" | wc -l)
n_ind=$(bcftools view -H "${PREFIX}_indels_raw.vcf.gz" | wc -l)
echo "records: total $n_total, snp $n_snp, indel $n_ind"
if [ $((n_snp + n_ind)) -ne "$n_total" ] && [ "${ALLOW_DROP:-0}" != 1 ]; then
    echo "ERROR: snp + indel ($((n_snp + n_ind))) != total ($n_total): MNP, mixed or reference-only records would be dropped or doubled." >&2
    echo "       Run: bcftools norm -m- --atomize <in> first, or set ALLOW_DROP=1 to accept the loss." >&2
    exit 2
fi

# The "|| INFO/X = \".\"" guard on every term replicates GATK's "missing => PASS": MQRankSum and ReadPosRankSum
# are undefined at hom-alt sites, so without it every hom-alt site silently fails and vanishes.
bcftools filter -i '
    QUAL >= 30 &&
    (INFO/QD >= 2.0 || INFO/QD = ".") &&
    (INFO/FS <= 60.0 || INFO/FS = ".") &&
    (INFO/MQ >= 40.0 || INFO/MQ = ".") &&
    (INFO/MQRankSum >= -12.5 || INFO/MQRankSum = ".") &&
    (INFO/ReadPosRankSum >= -8.0 || INFO/ReadPosRankSum = ".") &&
    (INFO/SOR <= 3.0 || INFO/SOR = ".")
' "${PREFIX}_snps_raw.vcf.gz" -Oz -o "${PREFIX}_snps_filtered.vcf.gz"

bcftools filter -i '
    QUAL >= 30 &&
    (INFO/QD >= 2.0 || INFO/QD = ".") &&
    (INFO/FS <= 200.0 || INFO/FS = ".") &&
    (INFO/ReadPosRankSum >= -20.0 || INFO/ReadPosRankSum = ".") &&
    (INFO/SOR <= 10.0 || INFO/SOR = ".")
' "${PREFIX}_indels_raw.vcf.gz" -Oz -o "${PREFIX}_indels_filtered.vcf.gz"

# SNP and indel records interleave by position, so plain concat aborts on multi-contig input; -a needs indexes.
bcftools index -f "${PREFIX}_snps_filtered.vcf.gz"
bcftools index -f "${PREFIX}_indels_filtered.vcf.gz"
bcftools concat -a "${PREFIX}_snps_filtered.vcf.gz" "${PREFIX}_indels_filtered.vcf.gz" | \
    bcftools sort -Oz -o "${PREFIX}_all_filtered.vcf.gz"
bcftools index -t -f "${PREFIX}_all_filtered.vcf.gz"
rm -f "${PREFIX}_snps_raw.vcf.gz" "${PREFIX}_indels_raw.vcf.gz"

n_snp_kept=$(bcftools view -H "${PREFIX}_snps_filtered.vcf.gz" | wc -l)
n_ind_kept=$(bcftools view -H "${PREFIX}_indels_filtered.vcf.gz" | wc -l)
n_out=$(bcftools view -H "${PREFIX}_all_filtered.vcf.gz" | wc -l)
echo "kept: snp $n_snp_kept of $n_snp, indel $n_ind_kept of $n_ind, total $n_out of $n_total"
if [ "$n_out" -eq 0 ]; then
    echo "ERROR: every record was filtered out; inspect the annotation distributions before relaxing thresholds." >&2
    exit 3
fi

tstv_before=$(bcftools stats "$IN" | grep '^TSTV' | cut -f5)
tstv_after=$(bcftools stats "${PREFIX}_all_filtered.vcf.gz" | grep '^TSTV' | cut -f5)
echo "Ti/Tv before: $tstv_before  after: $tstv_after  (expected about 2.0-2.1 for WGS, 3.0-3.3 for WES)"
case "${ASSAY:-}" in
    wgs) floor=2.0 ;;
    wes) floor=3.0 ;;
    *) floor= ;;
esac
if [ -n "$floor" ] && [ -n "$tstv_after" ] && awk -v t="$tstv_after" -v f="$floor" 'BEGIN{exit !(t < f)}'; then
    echo "WARNING: Ti/Tv $tstv_after is below $floor for ASSAY=${ASSAY}: filtering may be too loose (random errors have Ti/Tv about 0.5)." >&2
fi
echo "Output: ${PREFIX}_all_filtered.vcf.gz"
