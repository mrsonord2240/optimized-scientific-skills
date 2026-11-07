# SnpEff and ANNOVAR

No script; run the command and pin the engine version and transcript set.

```bash
# Human GRCh38 database needs 3-4 GB in memory: give the JVM 8 GB or more
snpEff -Xmx8g ann GRCh38.105 norm.vcf > out.vcf
snpEff -Xmx8g ann GRCh38.105 norm.vcf | SnpSift annotate clinvar.vcf.gz > annotated.vcf
```

```bash
table_annovar.pl norm.vcf humandb/ -buildver hg38 -out annotated -remove \
    -protocol refGene,gnomad30_genome,clinvar_20230416,dbnsfp42a \
    -operation g,f,f,f -nastring . -vcfinput
```

- SnpEff and SnpSift use a single dash for the version: `snpEff -version`.
- Normalize first (`bcftools norm -f ref.fa -m-any`). SnpEff `HIGH/MODERATE/LOW/MODIFIER` and ANNOVAR `exonic;splicing` are triage buckets, not evidence.
- The engines differ in splice-region width, up/downstream windows, HGVS 3'-shifting and severity ordering; do not compare a variant annotated by one with another. Check HGVS shifting per version.
- ANNOVAR protocols pair with operations: g gene, f filter, r region. The database versions in `-protocol` must match the build.
- Flatten SnpEff `ANN` with `routes/parse-csq.md`.

Done when the engine, version, database and transcript set are stated with the annotated VCF.
