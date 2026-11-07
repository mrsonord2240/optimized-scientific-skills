# Choosing a site-level method

No script: pick the row, then open the route it names.

| Dataset | Method | Route |
|---------|--------|-------|
| Single exome or gene panel, any non-model organism, or no truth resources | Hard filters, VETS or NVScoreVariants. Never VQSR | `routes/hard-filter.md` |
| One deep human WGS, or about 30 or more jointly genotyped exomes, with HapMap/Omni/Mills available | VQSR | `routes/vqsr.md` |
| Biobank or gnomAD scale | Allele-specific VQSR (`-AS`, `AS_*` annotations) or VETS | `routes/vqsr.md` for the plain form |
| Single sample where VQSR has too few variants | NVScoreVariants, then `FilterVariantTranches` | GATK docs for the installed build |
| Modern replacement for VQSR | VETS: `ExtractVariantAnnotations`, `TrainVariantAnnotationsModel`, `ScoreVariantAnnotations` (BETA; validate against a truth set first) | GATK docs for the installed build |
| DeepVariant or DRAGEN output | Filter on the caller's own QUAL and FILTER. Do not add GATK hard filters | none |
| Somatic tumor-normal | FilterMutectCalls plus caller-specific thresholds | `routes/somatic.md` |

- VQSR fits a Gaussian mixture over 6-8 annotations and needs tens of thousands of variants. On a single exome or panel it can report "converged" while VQSLOD is noise.
- On exomes, never give VQSR `-an DP`: capture depth tracks bait design, not truth.
- Confirm tool availability with `gatk --list` on the installed build before scripting VETS or NVScoreVariants. `CNNScoreVariants` is deprecated (replaced by NVScoreVariants); VETS is BETA.
- Rare-variant work: a real singleton looks like an artifact. Use tiered, looser thresholds with an orthogonal check (manual review, replication) instead of one blanket cutoff.

Done when you have told the user which method applies and why, and run its route.
