# bcftools annotation (rsIDs, frequencies, ClinVar, consequences)

```bash
bash scripts/annotate_vcf.sh input.vcf.gz output.vcf.gz --reference ref.fa --dbsnp dbsnp.vcf.gz \
    --gnomad gnomad.vcf.gz --gnomad-field fafmax_faf95_max --clinvar clinvar.vcf.gz --gff3 genes.gff3.gz
```

Every source is optional. Steps run in order: normalize (`--reference`), rsIDs, gnomAD value into a new `INFO/gnomAD_FREQ` tag, ClinVar `CLNSIG`/`CLNDN`, `bcftools csq -p a` (`--gff3`, needs `--reference`). Output is indexed and `output.summary.txt` gives the hit count of each step. The script refuses unindexed inputs, a source sharing no contig names with the input, and an INFO tag missing from a source header.

- The gnomAD value goes into a NEW tag. Copying into `INFO/AF` keeps the input's own AF wherever gnomAD has no record, so "absent from gnomAD" can no longer be detected. For v4 use the filtering AF `fafmax_faf95_max` (joint exome+genome: `fafmax_faf95_max_joint`), not global AF; see `routes/population-frequency.md`.
- `annotate -a` needs an indexed target, so each step writes an indexed file and cannot read a pipe.
- `csq -p r` (default) exits on the first unphased het. `-p a` takes GTs as is and merges nearby hets in a transcript into one consequence (a frameshift can absorb a downstream stop); `-p m` merges all GTs into one haplotype; `-p s` skips unphased hets, emitting no consequence. Use `-p a` for unphased short-read calls and say so.
- Other recipes (BED/TAB annotation, removing fields, `--set-id`, chromosome renaming, database downloads): `references/bcftools-recipes.md`.
- A ClinVar assertion is a research lead; classification is `variant-calling/clinical-interpretation`.

Done when the output VCF exists, the summary hit counts are read, and a zero-hit step is investigated (build or contig mismatch) rather than accepted.
