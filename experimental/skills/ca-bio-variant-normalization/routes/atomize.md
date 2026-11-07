# Atomize MNPs

```bash
bcftools norm --atomize input.vcf.gz -Oz -o atomized.vcf.gz
bcftools norm --atomize --old-rec-tag ORIGINAL input.vcf.gz -Oz -o atomized.vcf.gz   # keep the original record as INFO
```

`chr1 100 . ATG GCA` becomes `A G` at 100, `T C` at 101, `G A` at 102. Split multiallelics first (`routes/normalize.md` does it in order). Needs bcftools 1.12+; older versions use `vt decompose_blocksub`.

| Caller | Emits MNPs |
|--------|------------|
| FreeBayes | Yes |
| Octopus | Yes (block substitutions) |
| GATK HaplotypeCaller | No (may emit nearby SNPs in one haplotype block) |
| DeepVariant | Rarely |

Atomize to match across callers, dbSNP, ClinVar or gnomAD. Do not annotate consequence on the atomized copy: two adjacent SNVs in one codon can each look synonymous while the real MNV is missense or nonsense. Keep the atomized copy for matching and the un-atomized copy for annotation.

```bash
bcftools csq -p a -f ref.fa -g genes.gff3.gz in.vcf.gz
```

- `csq` defaults to `-p r` (phased only) and stops with `Unphased heterozygous genotype` on ordinary calls. For unphased input use `-p a` (takes GTs as given and builds haplotypes regardless of phase, merging nearby hets).
- `-p m` merges all GTs into one haplotype, trans-phased pairs included. `-p s` skips unphased hets, so it cannot create separate haplotypes; check the emitted BCSQ on your version.

Done when the atomized VCF exists and annotation was run on the un-atomized copy.
