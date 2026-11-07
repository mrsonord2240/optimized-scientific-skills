# Errors

| Error | Cause | Fix |
|-------|-------|-----|
| `REF does not match` | Wrong reference or build, contig naming | Use the FASTA the VCF was called against (`routes/ref-mismatch.md`) |
| `not sorted` | Unsorted input | `bcftools sort` first |
| Duplicate records after splitting | Same position twice | `bcftools norm -d exact` |
| `--atomize` unrecognized | bcftools below 1.12 | Upgrade, or `vt decompose_blocksub` |
| Spurious `ALT=*` records after atomize and split | Atomized before splitting | `norm -m- \| norm --atomize` |
| `Unphased heterozygous genotype` from `bcftools csq` | csq defaults to `-p r` | Set `-p a`, `-p m` or `-p s` (`routes/atomize.md`) |
| Split records carry the wrong per-allele AF or AD | Custom field declared `Number=.` | Declare `A`, `R` or `G` in the header |
| HGVS `c.` disagrees with VCF POS | 5' genomic left-align versus HGVS 3'-rule | Expected (`routes/hgvs.md`) |
| `bcftools index`: "not BGZF compressed" | A REF mismatch aborted the last pipe stage | Run `scripts/check_ref.sh`; use `pipefail` |

Done when the failing command has run to completion.
