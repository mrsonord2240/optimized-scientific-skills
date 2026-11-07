# REF allele mismatches

```bash
bash scripts/check_ref.sh reference.fa input.vcf.gz
```

Counts REF alleles that do not match the FASTA; exit 1 if any. Full list: `bcftools norm -f reference.fa -c w input.vcf.gz -Ou > /dev/null`.

`bcftools norm -c` modes: `e` error and exit (default), `w` warn and continue, `x` exclude mismatches, `s` set REF from the reference.

```bash
bcftools norm -f reference.fa -c x input.vcf.gz -Oz -o clean.vcf.gz   # drop mismatches
bcftools norm -f reference.fa -c s input.vcf.gz -Oz -o fixed.vcf.gz   # rewrite REF
```

- Mismatches usually mean a wrong build, `chr1` versus `1` contig naming, or a masked FASTA. Find the cause first; one off-by-one in a contig shifts every indel.
- `-c s` can hide a liftover error. When the cause is unknown, `-c x` is safer.

Done when `check_ref.sh` exits 0, or the mismatches are excluded or fixed with the cause named.
