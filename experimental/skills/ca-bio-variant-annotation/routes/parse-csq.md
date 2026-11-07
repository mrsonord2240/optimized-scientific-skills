# Flatten VEP CSQ or SnpEff ANN

```bash
python scripts/parse_csq.py annotated.vcf.gz --mane-only --impact HIGH,MODERATE --out parsed_annotations.tsv
```

One TSV row per transcript block: chrom, pos, ref, alt, then every field named in the VCF's CSQ format line (or the 15 standard SnpEff ANN fields). Needs `cyvcf2`. `--mane-only` keeps blocks with a populated `MANE_SELECT`; `--impact` keeps listed impact values. It refuses a VCF with neither CSQ nor ANN, and `--mane-only` where CSQ has no `MANE_SELECT` field.

- Prefer the MANE block over worst-consequence. Impact is a triage bucket.
- At the command line, `bcftools +split-vep -f '%CHROM\t%POS\t%SYMBOL\t%Consequence\n' -s worst out.vcf.gz` does the same; `-s` and `-p` choose which blocks surface.
- Fields come from the header, so a different VEP plugin set changes the columns.

Done when `parsed_annotations.tsv` exists and the number of variants and blocks is stated.
