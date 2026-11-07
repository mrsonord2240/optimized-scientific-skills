# Split and join multiallelic sites

```bash
bcftools norm -m-any input.vcf.gz -Oz -o split.vcf.gz      # split
bcftools norm -m+any input.vcf.gz -Oz -o merged.vcf.gz     # join biallelic records
```

`chr1 100 . A G,T ... GT 1/2` splits to `A G` with `1/0` and `A T` with `0/1`.

| Option | Effect |
|--------|--------|
| `-m-any` | Split all multiallelic sites |
| `-m-snps` / `-m-indels` | Split multiallelic SNPs / indels only |
| `-m-both` | Same as `-m-any` when splitting |
| `-m+any` | Join biallelic sites into multiallelic |
| `-m+snps` / `-m+indels` | Join SNPs / indels |
| `-m+both` | Join SNPs and indels into separate records (never merges a SNP with an indel) |

The type string (`any`/`both`) only matters when joining.

| Downstream tool | Split? |
|-----------------|--------|
| PLINK, PLINK2, most GWAS tools | Yes |
| ClinVar matching | Yes |
| Hail | No (handles multiallelics; splitting loses information) |
| `bcftools csq` | No |
| VEP | Either |

When a tool does not need splitting, keep multiallelic sites intact.

- Splitting loses that both alleles sit in one individual: 1/2 becomes 0/1 in both records. That breaks compound-het detection, and naive frequency calculation can double-count samples.
- Per-allele fields follow the header `Number`: `A` (AF, AC) takes the k-th value; `R` (AD) has REF first; `G` (PL, GL) needs the genotype-index formula; `Number=.` is carried whole onto every split record. A custom field mis-declared `.` when it is really `A` gives each record the wrong allele's value; fix the header.
- Join back (`-m+`) only if the split records were not independently filtered; removing one allele of a 1/2 site leaves a misleading record. Joining cannot always rebuild the original PLs.
- Spanning-deletion `*` alleles reference a deletion on another line; bcftools keeps them consistent, naive splitters corrupt them. Use `--keep-sum AD` to preserve summed allele depth across split records.
- Split once, early, and stay biallelic.

Done when no record has more than one ALT (or you have kept them deliberately) and per-allele fields were checked.
