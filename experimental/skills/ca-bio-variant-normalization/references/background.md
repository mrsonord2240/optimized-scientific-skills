# Background

Not needed to run any route.

## Normalization algorithm (as in `vt normalize` and `bcftools norm`)

- Right-trim: while all alleles end in the same base and all have length 2 or more, drop the last base of each.
- Left-extend: while any allele has length 0, prepend the reference base to the left of every allele and decrement POS.
- Left-trim: while all alleles start with the same base and all have length 2 or more, drop the first base and increment POS.

The right-trim then left-roll loop shifts an indel through a repeat to its leftmost equivalent position. Parsimony alone leaves an indel positionally ambiguous inside a repeat; left-alignment alone can leave redundant shared bases. Together they give one `(CHROM, POS, REF, ALT)` key, which is what makes join, merge and match safe.

## Why a one-base-off indel is a silent miss

In reference `...AAAAAAA...`, deleting any one of the seven A's gives the same alternate haplotype, so callers emit different POS. A non-left-aligned indel is a valid VCF line with a different tuple, so annotation keyed on the tuple misses the dbSNP, ClinVar or gnomAD entry at the canonical coordinate, and a clinically actionable variant is reported absent or novel with no error. `bcftools merge` of normalized and un-normalized cohorts emits duplicate rows and splits allele frequencies.

## Why split before atomize

`--atomize` defaults to `--atom-overlaps '*'`, so atomizing a multiallelic first rewrites a 1/2 SNV site as `A>C,*` and `A>G,*`, and the later split emits two spurious `A>*` records (checked on bcftools 1.21 and 1.24).
