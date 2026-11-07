---
name: ca-bio-variant-normalization
category: Data Analysis
description: Use when putting a VCF into canonical left-aligned, split, atomized form before comparing callers, merging cohorts or matching variants against databases.
tool_type: cli
primary_tool: bcftools
license: MIT
author: GPTomics
---

# Variant Normalization

Pick the one row that matches your task. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your task | Read |
|-----------|------|
| Normalize a VCF (split, atomize, left-align) for annotation, merging or matching | `routes/normalize.md` |
| Compare variants from two or more callers or cohorts | `routes/compare-callers.md` |
| Check REF alleles against the reference, or fix or drop mismatches | `routes/ref-mismatch.md` |
| Split or join multiallelic sites, per-allele fields, which tools need splitting | `routes/multiallelic.md` |
| Decompose MNPs and still annotate consequence correctly | `routes/atomize.md` |
| vt, GATK or bcftools disagree on counts | `routes/tool-discordance.md` |
| HGVS `c.` position disagrees with VCF POS | `routes/hgvs.md` |
| Estimate how much of a VCF needs normalizing | `routes/check.md` |
| Remove duplicate records, or prepare for PLINK | `routes/dedup-plink.md` |
| Annotate, filter, merge or call | `routes/other-skills.md` |
| An error message | `routes/errors.md` |

Rules on every route:

- Order is split (`-m-`), then atomize, then left-align with `-f ref.fa`. Atomizing a multiallelic first writes spurious `ALT=*` records.
- The FASTA must be the exact reference the VCF was called against and the one used downstream; a build, contig-name or masking mismatch gives wrong coordinates. Never paper over a REF mismatch with `-c w`.
- Use one tool and one flag set across every cohort you will compare, and record the command.
- An un-normalized indel is a valid VCF line that silently misses its dbSNP, ClinVar or gnomAD entry; no error is thrown.
- Normalization is idempotent: normalizing a normalized VCF changes nothing.

Tested with bcftools 1.24 (`--atomize` needs 1.12+), cyvcf2 0.30+.
