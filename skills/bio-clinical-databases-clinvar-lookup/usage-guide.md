# ClinVar Lookup - Usage Guide

## Overview

Query ClinVar for variant pathogenicity classifications, ClinGen Variant Curation Expert Panel (VCEP) curations, and germline-vs-somatic interpretations against the ACMG/AMP framework with ClinGen SVI specifications. The skill addresses the VCV/SCV/RCV identifier hierarchy, the 2024 XML schema overhaul, star-rating override semantics, and the ClinGen Allele Registry as the canonical cross-database join key.

## Prerequisites

```bash
pip install requests cyvcf2 pandas lxml

# Local VCF (weekly snapshot; pin to monthly archive for reproducibility)
mkdir -p clinvar/$(date +%Y%m); cd clinvar/$(date +%Y%m)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# For reproducibility, use the monthly archive
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2025/clinvar_20251101.vcf.gz
```

## Quick Start

Tell the agent what to do:
- "Look up the ClinVar classification for BRAF V600E and report star rating and VCEP curation status"
- "Annotate my somatic VCF with germline + oncogenicity + somatic-clinical-impact classifications from ClinVar 2024 schema"
- "Resolve these variants to ClinGen Allele Registry CA IDs for cross-database joining with gnomAD"
- "Triangulate conflicting interpretations for this VUS: list all SCVs with submitter and star, identify the highest-confidence assertion"
- "Filter my list to retain only variants with star >=2 OR ClinGen VCEP curation"

## Example Prompts

### Single Variant Queries

> "What is the ClinVar VCV classification for chr17:43106487:A:C in BRCA1, and is there an ENIGMA-VCEP curation? Report germline classification, last evaluated date, and review status."

> "Look up rs121913529 in ClinVar; show all SCV-level assertions with submitter, classification, condition, and date last evaluated."

> "Resolve chr7:140453136:A:T to ClinGen Allele Registry CA ID and report the ClinVar VariationID linkage."

### Condition-Stratified Queries

> "For BRCA2 c.5946delT, give me the RCV-level classification for each condition rather than the VCV aggregate."

> "Find all ClinVar variants in TP53 with germline Pathogenic + somatic Tier I oncogenicity classifications under the 2024 tripartite schema."

### Conflict Resolution and Star Filtering

> "List all conflicting-interpretation variants in BRCA1, group by submitter pair, and identify whether the conflict is P/LP vs LB/B (clinically meaningful) or P vs LP (often immaterial)."

> "Filter my variant list to keep only star >=2 OR ClinGen VCEP-curated; flag stale assertions older than 36 months."

### Bulk Annotation

> "Annotate this exome VCF with ClinVar germline + oncogenicity + somatic-clinical-impact + review status + VCEP affiliation using bcftools."

> "For my 50k rare-variant cohort, batch-query ClinVar via myvariant.info with `fields=clinvar.rcv.review_status,clinvar.variant_id` and merge against gnomAD grpmax FAF95."

### Cross-Database Join

> "Resolve these 200 HGVS-g variants to CA IDs and then look each up in ClinVar, gnomAD v4 exomes, COSMIC, and MAVEdb in parallel."

## What the Agent Will Do

1. Parse the input identifier (rsID / HGVS / chromosomal coords / gene-protein) and resolve to a canonical form (CA ID or VCV).
2. Choose query mode by batch size: REST esummary for <10; myvariant.info for 10-1000; local VCF or bulk XML for >1000.
3. Pull VCV-level aggregate AND condition-stratified RCV-level classifications when the user references a specific phenotype.
4. Read the 2024 tripartite (Germline / SomaticClinicalImpact / Oncogenicity) classifications separately.
5. Apply the star-rating override hierarchy (4 > 3 > 2 > 1 > 0); flag stale assertions and conflicting interpretations.
6. Cross-check pathogenic calls against gnomAD grpmax FAF95 for BS1/BA1 reconciliation.
7. Return CA ID + VCV + RCV + per-classification labels + review status + last evaluated date.

## Tips

- Use the monthly archive (first Thursday of each month) for reproducible analyses; weekly releases are not archived.
- `CLNSIG` in `clinvar.vcf.gz` is VCV-level (variant-level aggregate); for condition-specific classification parse RCV-level XML.
- 2024 XML schema replaces `<ClinVarSet>` with `<VariationArchive>`; pipelines built before September 2024 must be re-targeted.
- Star rating 3 (VCEP) supersedes lower-star records per ClinGen FDA Recognition 2018; do not auto-aggregate by date.
- `criteria provided, conflicting classifications` (formerly "conflicting interpretations") is 1-star (not 2-star as sometimes reported); inspect `CLNSIGCONF` to see whether the conflict is clinically meaningful.
- Use ClinGen Allele Registry CA ID (`https://reg.clinicalgenome.org/`) for any cross-database join; ClinVar VariationID was renumbered during the 2017 schema redesign.
- ClinVar somatic classifications (`ONCDN`, `SCIDN`) were added in 2024; pre-2024 pipelines miss them silently.
- Conflict resolution is slow: only ~4% of BRCA1 missense VUS conflicts have reached consensus despite years of effort.
- For pathogenicity classification logic (PVS1 decision tree, Pejaver 2022 calibrated PP3/BP4 thresholds, Tavtigian point system), defer to `clinical-databases/acmg-classification`; this skill is for querying ClinVar, not classification.

## Reconciliation: When Sources Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| ClinVar P vs gnomAD AF > 1% | Variant is true founder allele in unstratified gnomAD subset, OR ClinVar P is a stale low-star assertion | Check `grpmax_faf95` excluding bottleneck groups; check ClinVar star rating |
| ClinVar P vs AlphaMissense < 0.1 | Variant in NMD-escape region, alternative isoform, or ClinVar P is mis-curated | Check Pejaver 2022 calibration in `acmg-classification` skill; cross-check VCEP |
| VCEP 3-star P vs commercial-lab 1-star B | VCEP review takes precedence | Weight the VCEP assertion; flag the submitter record as discordant |
| ClinVar VCV-level P vs RCV-level VUS for actual condition | VCV averages across conditions | Report at RCV level for condition-specific work |
| ClinVar P vs LOVD/HGMD discordant | LOVD/HGMD use different classification systems; HGMD "DM" != ACMG P | Triangulate against published evidence; do not auto-translate labels |
| ClinVar P missing for a known disease variant | Submission lag (~6-12 months typical for new findings) | Check published literature; flag for ClinVar submission |

## ClinVar Somatic vs Germline: 2024 Tripartite

The 2024 schema separates three orthogonal classifications, each with its own `ReviewStatus` and `DateLastEvaluated`:

- **GermlineClassification:** Pathogenic / Likely Pathogenic / VUS / LB / B per ACMG/AMP 2015 + SVI.
- **SomaticClinicalImpact:** Tier I / II / III / IV per AMP/ASCO/CAP 2017 (Li 2017 *J Mol Diagn*).
- **OncogenicityClassification:** Oncogenic / Likely Oncogenic / VUS / Likely Benign / Benign per ClinGen/CGC/VICC 2022 oncogenicity framework.

A single VCV can carry all three with distinct evaluations; the legacy "Pathogenic" label is now ambiguous if not qualified by classification type.

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why is this pathogenic variant 1-star?" | We report the star rating per record; stars describe review depth, and ClinVar entries are research evidence, not a clinical classification. |
| "ClinVar says P but gnomAD AF = 2%" | Reconciled via Whiffin FAF95 max-credible-AF framework; bottleneck-group rule applied. |
| "This VCV count differs from ClinVar.gov" | We pulled from the monthly archive (first-Thursday-of-month) for reproducibility; the live web is post-most-recent-weekly. |
| "Why wasn't the somatic variant flagged?" | Pre-2024 XML schema had no separate somatic field; we now read `ONCDN`/`SCIDN`/`SomaticClinicalImpact` per v2 schema. |
| "VarSome says LP but this says VUS" | Tool-specific aggregation rule differences; VarSome auto-applies PP3+PM2 by default per Tavtigian point system; we apply VCEP-specific PP3 calibration per CSpec. |
| "rsID match returned wrong variant" | rsID is a cluster identifier; multi-allelic rsIDs require allele-level resolution; we use SPDI or CA ID. |
| "Why retest a 2022-curated variant?" | Classifications drift as evidence accrues; ClinGen recommends annual re-review for active diagnostic variants. |

## Related Skills

- clinical-databases/acmg-classification - ACMG/AMP framework, Pejaver PP3/BP4 calibration, PVS1 decision tree
- clinical-databases/myvariant-queries - Multi-database aggregation with ClinVar overlay
- clinical-databases/variant-prioritization - Rare-disease filtering pipeline
- clinical-databases/gnomad-frequencies - Population frequency for BS1/BA1 cross-check
- variant-calling/clinical-interpretation - Clinical reporting workflow
