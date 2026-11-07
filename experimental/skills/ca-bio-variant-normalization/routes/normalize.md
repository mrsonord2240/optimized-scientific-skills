# Normalize a VCF

```bash
bash scripts/normalize_vcf.sh reference.fa input.vcf.gz normalized.vcf.gz
ATOMIZE=0 bash scripts/normalize_vcf.sh reference.fa input.vcf.gz normalized.vcf.gz   # no MNPs in the input
```

The script checks REF alleles first and stops with a clear message on a mismatch, then runs `bcftools norm -m-`, `norm --atomize`, `norm -f ref.fa` in that order, indexes, and prints record counts before, after split plus left-align, and after atomizing. `ATOMIZE=0` skips atomizing for callers that emit no MNPs (GATK HaplotypeCaller).

Normalize before: comparing callers, annotating against dbSNP, ClinVar or gnomAD, `bcftools merge` (otherwise duplicate rows and split allele frequencies), and `bcftools isec`. It is safe to skip only when one caller made every variant and nothing is compared or looked up.

- A variant is normalized when it is parsimonious (no base trimmable from either end; indels keep one anchor base, so a deletion of `A` is `REF=CA ALT=C`) and left-aligned (POS cannot move left).
- Homopolymer example, reference `AAAAAA` at chr1:100: `100 AA A` is canonical; `101 AA A` and `104 AA A` are the same deletion written further right and will not match.
- Pipe form: `bcftools norm -m- in.vcf.gz | bcftools norm --atomize | bcftools norm -f ref.fa -Oz -o out.vcf.gz`; use `set -o pipefail`, since a REF mismatch aborts the last stage after writing a 0-record, non-BGZF file and `bcftools index` then fails with "not BGZF compressed".
- Single pass `bcftools norm -f ref.fa -m-any --atomize` also splits before atomizing. `-m-any` alone skips MNP atomization.
- Reference must be indexed (`samtools faidx reference.fa`).

Done when `normalized.vcf.gz` and its index exist and the script printed its counts.
