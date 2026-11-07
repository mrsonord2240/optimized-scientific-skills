# Allele harmonisation

## Allele Harmonisation (Critical Pre-Step)

Mismatched effect alleles silently invert signs of betas, collapsing PP.H4 into PP.H3. Required steps before coloc:

1. Merge GWAS and eQTL summary stats by SNP ID (rsID or chr:pos:ref:alt).
2. Mark SNP-pairs as `same` (A1/A2 match) or `flip` (A1/A2 swap); for non-palindromic pairs reported on opposite strands (e.g. A/G vs T/C), complement dataset 2's alleles first, then classify; drop SNPs that match neither.
3. For `flip` rows, negate the second dataset's beta (and swap A1/A2).
4. Drop palindromic SNPs (A/T or C/G) at MAF > 0.42; their strand cannot be inferred from coding alone (TwoSampleMR `harmonise_data` standard cutoff).
5. Verify genome build alignment (hg19 vs hg38 must match; lift over if not).

```bash
Rscript scripts/harmonise.R gwas.tsv eqtl.tsv harmonised.tsv     # columns SNP, A1, A2, BETA, MAF
# or in R: source('scripts/harmonise.R'); m <- harmonise(df1, df2)
```

`harmonise()` (`scripts/harmonise.R`) merges on SNP, resolves same / flip / strand-complement coding for non-palindromic pairs, negates `BETA.2` (and sets `MAF.2 = 1 - MAF.2`) on flips, and drops unresolvable pairs and palindromic SNPs at MAF > 0.42. Output columns are suffixed `.1` / `.2`.

Harmonisation pitfalls to watch for:

- **Allele coding mismatch.** GWAS may report effect allele as A1 while eQTL reports it as A2. Always check both and flip betas where needed.
- **Build mismatch.** hg19 GWAS coords + hg38 eQTL coords silently merge on rsID but break on chr:pos. Lift over with `rtracklayer::liftOver` or CrossMap before merging.
- **Strand mismatch.** Non-palindromic SNPs coded on opposite strands (A/G vs T/C) are resolved by the complement step in `harmonise()`. Palindromic SNPs cannot be resolved this way (next bullet).
- **Palindromic SNPs at high MAF.** A/T and C/G SNPs at MAF > 0.42 cannot be unambiguously strand-resolved; drop them or resolve with reference-panel MAF.
- **Multi-allelic SNPs.** Many summary stats collapse multi-allelic loci by keeping only the most-frequent alt; if datasets pick different alts, harmonisation drops the SNP. Split on chr:pos:ref:alt as a unique key.
- **rsID dependence.** rsID can be remapped across dbSNP builds (e.g. merge of two rsIDs into one). Prefer chr:pos:ref:alt keys for cross-study merges.
