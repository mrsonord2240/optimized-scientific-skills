# bcftools annotate recipes

Read when `routes/annotate-vcf.md` is not enough. `annotate -a <vcf>` needs both the target and the source indexed; a chr or build mismatch exits 0 with no matches, so compare `bcftools index -s` of both first.

```bash
# Copy fields from a VCF source
bcftools annotate -a dbsnp.vcf.gz -c ID input.vcf.gz -Oz -o with_rsids.vcf.gz
bcftools annotate -a gnomad.vcf.gz -c INFO/gnomAD_FREQ:=INFO/fafmax_faf95_max in.vcf.gz -Oz -o out.vcf.gz   # new tag

# BED source (4 columns: CHROM, START, END, NAME) and TAB source (CHROM, POS, VALUE): the header line is required
echo '##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region">' > header.txt
bcftools annotate -a regions.bed.gz -c CHROM,FROM,TO,INFO/REGION -h header.txt in.vcf.gz -Oz -o out.vcf.gz
bcftools annotate -a scores.tab.gz -c CHROM,POS,INFO/SCORE -h header.txt in.vcf.gz -Oz -o out.vcf.gz

# Remove fields
bcftools annotate -x INFO/DP,INFO/MQ in.vcf.gz -Oz -o clean.vcf.gz
bcftools annotate -x FORMAT/AD,FORMAT/PL in.vcf.gz -Oz -o clean.vcf.gz
bcftools annotate -x ^INFO/DP,INFO/AF in.vcf.gz -Oz -o minimal.vcf.gz      # keep only these INFO fields

# Variant IDs; the + prefix sets only missing IDs
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' in.vcf.gz -Oz -o with_ids.vcf.gz
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' in.vcf.gz -Oz -o with_ids.vcf.gz

# Rename contigs: rename.txt is old<TAB>new (chr1 1; chrX X; chrM MT, or the reverse)
bcftools annotate --rename-chrs rename.txt in.vcf.gz -Oz -o renamed.vcf.gz

# Restrict to regions for large sources
bcftools annotate -a database.vcf.gz -c INFO/AF -r chr1,chr2 in.vcf.gz -Oz -o out.vcf.gz
```

Column specifiers: `ID`, `QUAL`, `FILTER`, `INFO`, `INFO/TAG`, `+INFO/TAG` (append), `-INFO/TAG` (remove).

Troubleshooting: "no such tag" means the source header lacks it (`bcftools view -h src.vcf.gz | grep '^##INFO'`); "index required" means `tabix -p vcf src.vcf.gz`; "chromosome not found" means a naming mismatch (rename as above).

Databases: dbSNP `GCF_000001405.40.gz` (ftp.ncbi.nih.gov/snp/latest_release/VCF); ClinVar `clinvar.vcf.gz` (ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38) with `CLNSIG`, `CLNDN`; gnomAD v4 per-chromosome `gnomad.genomes.v4.0.sites.chrN.vcf.bgz` (fields `AF`, `AF_grpmax`, `fafmax_faf95_max`; v4 renamed popmax to grpmax). Download only the chromosomes needed.

bcftools csq output: `BCSQ=consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change`. Query with `bcftools query -f '%CHROM\t%POS\t%INFO/BCSQ\n'`; filter with `bcftools view -i 'INFO/BCSQ~"stop_gained" || INFO/BCSQ~"frameshift"'`.
