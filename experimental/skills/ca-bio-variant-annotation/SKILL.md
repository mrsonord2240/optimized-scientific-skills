---
name: ca-bio-variant-annotation
category: Data Analysis
description: Use when annotating a human VCF with consequences, population frequencies and pathogenicity evidence, or when choosing the transcript set, engine and predictors for that annotation.
tool_type: mixed
primary_tool: VEP
license: MIT
author: GPTomics
---

# Variant Annotation

Pick the one row that matches your task. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your task | Read |
|-----------|------|
| Annotate on MANE Select with Ensembl VEP (clinical default, GRCh38) | `routes/annotate-vep.md` |
| rsIDs, gnomAD frequency, ClinVar or bcftools csq consequences with bcftools only | `routes/annotate-vcf.md` |
| SnpEff or ANNOVAR | `routes/snpeff-annovar.md` |
| Flatten VEP CSQ or SnpEff ANN into a table | `routes/parse-csq.md` |
| Which transcript to report on; why `--pick` is dangerous; HGVS and left-alignment | `routes/transcript-choice.md` |
| Does a stop, frameshift or splice consequence support PVS1 (NMD rule) | `routes/consequence-nmd.md` |
| Which predictor to use for PP3/BP4 (REVEL, AlphaMissense, CADD, SpliceAI) | `routes/predictors.md` |
| What allele frequency to use; gnomAD version; is "in gnomAD" benign | `routes/population-frequency.md` |
| Left-align, split multiallelics, final ACMG classification, filtering | `variant-calling/variant-normalization`, `variant-calling/clinical-interpretation`, `bio-variant-calling-filtering-best-practices` |

Rules on every route:

- A consequence belongs to the tuple (variant, transcript model, engine, engine version, parameters). On any discordant result, ask which transcript and tool produced it. Pin and record genome build, transcript set, engine and version, predictor and version, gnomAD version.
- Normalize before annotating (`bcftools norm -f ref.fa -m-any`); the scripts do. Never hand-derive HGVS from POS.
- Confirm the genome build and the chromosome naming of the VCF and every annotation source match. gnomAD v2 is GRCh37; v3 and v4 are GRCh38; MANE is GRCh38 only. A chr1 versus 1 mismatch annotates nothing and exits 0.
- Report on MANE Select (plus MANE Plus Clinical where assigned), not worst-consequence across transcripts.
- Impact buckets (HIGH, MODERATE) are triage, not ACMG evidence.
- This Skill supplies consequence, NMD, predictor and frequency inputs; it does not combine them into a classification.

Tested with bcftools 1.19+ and VEP 110+ command forms (not rerun in this recut), SnpEff 5.2+, ANNOVAR 2020Jun07+.
