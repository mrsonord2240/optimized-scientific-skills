# vt, bcftools and GATK disagree

| Behavior | `vt` | `bcftools norm` | GATK `LeftAlignAndTrimVariants` |
|----------|------|-----------------|--------------------------------|
| Left-align and parsimony (biallelic) | `normalize` | `-f` | yes |
| Split multiallelic | `vt decompose` | `-m-any` | `--split-multi-allelics` |
| MNP to SNPs | `vt decompose_blocksub` does it by default | only with `--atomize` (1.12+) | no |
| Block substitution or complex | `decompose_blocksub` | `--atomize` | limited |

They agree on simple biallelic indels and differ on decomposition: bcftools keeps an MNP as one record, vt emits separate SNPs, so the variant counts differ. If cohort A is atomized with vt and cohort B is left alone by bcftools, every MNP becomes a cohort-private variant.

- Use one tool and one exact flag set for every cohort you compare; record the command. `bcftools norm` is the production standard (split, atomize and left-align in one tool).
- GATK left-aligns only within `--max-leading-bases`; indels inside longer repeats may not fully left-align. Check `gatk LeftAlignAndTrimVariants --help` for the installed default before trusting STR-embedded indels.

Done when both cohorts used the same tool and flags, with the command recorded.
