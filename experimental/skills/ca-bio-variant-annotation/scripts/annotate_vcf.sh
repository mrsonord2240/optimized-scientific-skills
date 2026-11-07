#!/bin/bash
# Purpose: lightweight bcftools annotation of a VCF: normalize, add rsIDs, gnomAD frequency into a NEW tag,
#          ClinVar assertions, and bcftools csq consequences. Each source is optional.
# Usage:   scripts/annotate_vcf.sh input.vcf.gz output.vcf.gz [--reference ref.fa] [--dbsnp dbsnp.vcf.gz]
#            [--gnomad gnomad.vcf.gz] [--gnomad-field fafmax_faf95_max] [--clinvar clinvar.vcf.gz]
#            [--gff3 genes.gff3.gz]
#          --reference normalizes first (left-align, split multiallelics); --gff3 needs --reference and runs bcftools csq -p a.
# Output:  output.vcf.gz (+ .csi) and output.summary.txt with the hit count of every step.
# Refuses: missing or unindexed inputs, a source sharing no contig names with the input (chr1 vs 1 would annotate
#          nothing and exit 0), a requested INFO tag absent from the source header.
set -euo pipefail

die() { echo "REFUSED: $*" >&2; exit 2; }
command -v bcftools >/dev/null 2>&1 || die "bcftools not found on PATH"

[ $# -ge 2 ] || die "usage: $0 input.vcf.gz output.vcf.gz [--reference ref.fa] [--dbsnp f] [--gnomad f] [--gnomad-field TAG] [--clinvar f] [--gff3 f]"
INPUT=$1; OUTPUT=$2; shift 2
REFERENCE=""; DBSNP=""; GNOMAD=""; GNOMAD_FIELD="fafmax_faf95_max"; CLINVAR=""; GFF3=""
while [ $# -gt 0 ]; do
  case "$1" in
    --reference) REFERENCE=$2; shift 2;;
    --dbsnp) DBSNP=$2; shift 2;;
    --gnomad) GNOMAD=$2; shift 2;;
    --gnomad-field) GNOMAD_FIELD=$2; shift 2;;
    --clinvar) CLINVAR=$2; shift 2;;
    --gff3) GFF3=$2; shift 2;;
    *) die "unknown option $1";;
  esac
done
[ -z "$GFF3" ] || [ -n "$REFERENCE" ] || die "--gff3 needs --reference"

for f in "$INPUT" "$DBSNP" "$GNOMAD" "$CLINVAR"; do
  [ -z "$f" ] && continue
  [ -f "$f" ] || die "file not found: $f"
  bcftools index -s "$f" >/dev/null 2>&1 || die "indexed VCF required: $f (run: bcftools index -f \"$f\")"
done
[ -z "$REFERENCE" ] || [ -f "$REFERENCE" ] || die "reference not found: $REFERENCE"
[ -z "$GFF3" ] || [ -f "$GFF3" ] || die "gff3 not found: $GFF3"

check_contigs() {  # annotate exits 0 when no contig names match
  local shared
  shared=$(comm -12 <(bcftools index -s "$INPUT" | cut -f1 | LC_ALL=C sort -u) \
                    <(bcftools index -s "$1" | cut -f1 | LC_ALL=C sort -u))
  [ -n "$shared" ] || die "$1 shares no contig names with the input (chr prefix or build mismatch)"
}
check_tag() {  # $1 file, $2 INFO tag
  bcftools view -h "$1" | grep -q "^##INFO=<ID=$2," || die "INFO tag $2 not in the header of $1"
}

SUMMARY="${OUTPUT%.vcf.gz}.summary.txt"
: > "$SUMMARY"
WORK=$(mktemp -d); trap 'rm -rf "$WORK"' EXIT
CUR="$INPUT"
step=0
next() { step=$((step + 1)); echo "$WORK/step$step.vcf.gz"; }

if [ -n "$REFERENCE" ]; then
  NXT=$(next)
  bcftools norm -f "$REFERENCE" -m-any "$CUR" -Oz -o "$NXT"; bcftools index -f "$NXT"; CUR=$NXT
  echo "normalized (left-aligned, multiallelics split): $(bcftools view -H "$CUR" | wc -l | tr -d ' ') records" >> "$SUMMARY"
elif bcftools view -H "$CUR" | awk -F'\t' '$5 ~ /,/ {f=1; exit} END {exit !f}'; then
  echo "WARNING: multiallelic records present and no --reference given, so they were not split or left-aligned; annotate after normalizing." | tee -a "$SUMMARY" >&2
fi

if [ -n "$DBSNP" ]; then
  check_contigs "$DBSNP"
  NXT=$(next)
  bcftools annotate -a "$DBSNP" -c ID "$CUR" -Oz -o "$NXT"; bcftools index -f "$NXT"; CUR=$NXT
  echo "rsIDs from dbSNP: $(bcftools view -H "$CUR" | awk -F'\t' '$3 != "."' | wc -l | tr -d ' ') of $(bcftools view -H "$CUR" | wc -l | tr -d ' ') have an ID" >> "$SUMMARY"
fi

if [ -n "$GNOMAD" ]; then
  check_contigs "$GNOMAD"; check_tag "$GNOMAD" "$GNOMAD_FIELD"
  NXT=$(next)
  # A NEW tag: copying into INFO/AF keeps the input's own AF where gnomAD has no record, so "absent" is lost.
  bcftools annotate -a "$GNOMAD" -c "INFO/gnomAD_FREQ:=INFO/$GNOMAD_FIELD" "$CUR" -Oz -o "$NXT"; bcftools index -f "$NXT"; CUR=$NXT
  hits=$(bcftools query -f '%INFO/gnomAD_FREQ\n' "$CUR" | awk '$1 != "."' | wc -l | tr -d ' ')
  echo "gnomAD $GNOMAD_FIELD -> INFO/gnomAD_FREQ: $hits of $(bcftools view -H "$CUR" | wc -l | tr -d ' ') have a value; the rest are absent or not callable there, not necessarily rare" >> "$SUMMARY"
fi

if [ -n "$CLINVAR" ]; then
  check_contigs "$CLINVAR"; check_tag "$CLINVAR" CLNSIG
  NXT=$(next)
  bcftools annotate -a "$CLINVAR" -c INFO/CLNSIG,INFO/CLNDN "$CUR" -Oz -o "$NXT"; bcftools index -f "$NXT"; CUR=$NXT
  echo "ClinVar CLNSIG: $(bcftools query -f '%INFO/CLNSIG\n' "$CUR" | awk '$1 != "."' | wc -l | tr -d ' ') records have an assertion (research lead, not a classification)" >> "$SUMMARY"
fi

if [ -n "$GFF3" ]; then
  NXT=$(next)
  # -p a: take GTs as is; the default -p r exits on the first unphased het. Nearby hets in one transcript merge into one consequence.
  bcftools csq -p a -f "$REFERENCE" -g "$GFF3" "$CUR" -Oz -o "$NXT"; bcftools index -f "$NXT"; CUR=$NXT
  echo "bcftools csq: $(bcftools query -f '%INFO/BCSQ\n' "$CUR" | awk '$1 != "."' | wc -l | tr -d ' ') records have a BCSQ consequence" >> "$SUMMARY"
fi

cp "$CUR" "$OUTPUT"
bcftools index -f "$OUTPUT"
cat "$SUMMARY"
echo "wrote $OUTPUT and $SUMMARY"
