# Aligned DNA logo

```bash
Rscript scripts/logo_fasta.R aligned.fa logo.pdf bits
```

Input must contain one complete sequence per FASTA line, equal length, with A/C/G/T only; the source reader does not join wrapped FASTA records. This route uses ggseqlogo's sequence-vector bits mode and displays N. For explicitly supplied background use `Rscript scripts/logo_fasta.R aligned.fa logo.pdf relative 0.25 0.25 0.25 0.25` with positive A,C,G,T frequencies; the helper normalizes them and uses unsmoothed maximum-likelihood probabilities, not ggseqlogo's small-sample correction. The output can overwrite an existing file. Treat N below five as exploratory and disclose N and correction for all small samples; the source suggests 20 as a practical minimum, not a guarantee.

Done when alignment, alphabet, encoding, N, background, correction and rendered figure are checked.
