#!/bin/bash
# Reference: bcftools 1.19+ | Verify API if version differs
# Annotate VCF with rsIDs and population frequencies

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 <input.vcf.gz> <dbsnp.vcf.gz> <output.vcf.gz>"
    echo ""
    echo "Optional: Set GNOMAD_VCF environment variable to add population frequencies"
    exit 1
fi

INPUT="$1"
DBSNP="$2"
OUTPUT="$3"

# `annotate -a` needs both the target and source VCFs to be indexed. Check
# this up front: otherwise a missing index is reported later and can look like
# an annotation failure.
if ! command -v bcftools >/dev/null 2>&1; then
    echo "Error: bcftools is required but was not found on PATH" >&2
    exit 127
fi

# Check inputs and their indexes.
for vcf in "$INPUT" "$DBSNP"; do
    if [ ! -f "$vcf" ]; then
        echo "Error: File not found: $vcf" >&2
        exit 1
    fi
    if ! bcftools index -s "$vcf" >/dev/null 2>&1; then
        echo "Error: Indexed VCF required: $vcf (run: bcftools index -f \"$vcf\")" >&2
        exit 1
    fi
done

echo "Adding rsIDs from dbSNP..."
if [ -n "${GNOMAD_VCF:-}" ]; then
    if [ ! -f "$GNOMAD_VCF" ]; then
        echo "Error: GNOMAD_VCF was set but no file exists at: $GNOMAD_VCF" >&2
        exit 1
    fi
    if ! bcftools index -s "$GNOMAD_VCF" >/dev/null 2>&1; then
        echo "Error: Indexed gnomAD VCF required: $GNOMAD_VCF (run: bcftools index -f \"$GNOMAD_VCF\")" >&2
        exit 1
    fi
    # `annotate` exits 0 when every contig name differs (for example chr1 vs
    # 1), so fail before producing an apparently successful zero-hit result.
    shared_contigs=$(comm -12 \
        <(bcftools index -s "$INPUT" | cut -f1 | LC_ALL=C sort -u) \
        <(bcftools index -s "$GNOMAD_VCF" | cut -f1 | LC_ALL=C sort -u))
    if [ -z "$shared_contigs" ]; then
        echo "Error: Input and gnomAD VCFs have no shared contig names; check build and chr-prefix conventions." >&2
        exit 1
    fi
    echo "Adding population frequencies from gnomAD (INFO/gnomAD_AF)..."
    # annotate -a <vcf> needs an indexed target, so it cannot read the previous step from a pipe.
    IDS="${OUTPUT%.vcf.gz}.ids.vcf.gz"
    bcftools annotate -a "$DBSNP" -c ID "$INPUT" -Oz -o "$IDS"
    bcftools index -f "$IDS"
    # New tag: copying into INFO/AF would keep the input's own AF where gnomAD has no record.
    bcftools annotate -a "$GNOMAD_VCF" -c INFO/gnomAD_AF:=INFO/AF "$IDS" -Oz -o "$OUTPUT"
    rm -f "$IDS" "$IDS.csi"
else
    bcftools annotate -a "$DBSNP" -c ID "$INPUT" -Oz -o "$OUTPUT"
fi

bcftools index "$OUTPUT"

# Report
TOTAL=$(bcftools view -H "$OUTPUT" | wc -l | tr -d ' ')
WITH_ID=$(bcftools view -H "$OUTPUT" | awk -F'\t' '$3!="."' | wc -l | tr -d ' ')
PERCENT_WITH_ID=$(awk -v with_id="$WITH_ID" -v total="$TOTAL" 'BEGIN {
    if (total == 0) {
        print "0.0"
    } else {
        printf "%.1f", 100 * with_id / total
    }
}')

echo ""
echo "=== Annotation Complete ==="
echo "Output: $OUTPUT"
echo "Total variants: $TOTAL"
echo "With rsID: $WITH_ID (${PERCENT_WITH_ID}%)"
