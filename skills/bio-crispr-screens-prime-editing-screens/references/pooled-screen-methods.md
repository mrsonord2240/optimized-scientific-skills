# Pooled-Screen Methods: PRIME, MOSAIC and PE-BE Cross-Validation

Moved from SKILL.md. Read when planning a pooled PE screen, a saturation-mutagenesis library, or a BE-PE concordance analysis.

## PRIME Pooled Screen Methodology

**Ren X et al 2023 *Mol Cell* 83:4633** established the PRIME pooled prime-editing screen methodology (earlier 2023 bioRxiv preprint):

- pegRNA library covering thousands of intended variants, screened for specificity at design time (Ren 2023 used GuideScan2; add PRIDICT2 efficiency prediction for new designs)
- Lentiviral delivery in a PE-expressing cell line (Ren 2023: MOI 0.3 for the MYC-enhancer screen, MOI 0.5 for the variant screens)
- Selection on integration marker
- Time-course screen for variant function (e.g., drug sensitivity)
- Endpoint amplicon sequencing of each pegRNA target locus
- CRISPResso2 quantification of intended-edit %
- MAGeCK / drugZ-style hit calling on edit-efficient pegRNAs

**Quantified scale:** ~3,699 ClinVar variants installed in a single PRIME screen, alongside 1,304 breast-cancer GWAS variants.

## MOSAIC In Situ Saturation Mutagenesis

**MOSAIC (Hsu 2024, bioRxiv)** is a high-throughput in-situ saturation-mutagenesis prime-editing method with multiplexed read-out:

- Tile pegRNAs across protein domains for systematic mutagenesis
- Saturation: every possible amino acid change in a region
- Identify drug-resistance variants in real-time
- Smaller per-variant cell numbers (more variants total)

**Use case:** Cancer-drug-resistance variant scanning; protein-domain function mapping.

## Cross-Validate PE with Base Editor Screens

**Goal:** Confirm variant-function calls from PE with orthogonal BE screens.

**Approach:** Design parallel BE library for the same variants; run both screens; intersect hits.

```bash
# Inner-joins on variant_id; high_confidence = both FDR < 0.05 and the same sign of LFC.
# Columns needed: be_fdr, be_lfc (BE file) and pe_fdr, pe_lfc (PE file).
python scripts/crossvalidate_pe_be.py be_screen_hits.tsv pe_screen_hits.tsv --fdr 0.05
```

**Critical:** PE-only hits in BE-coverable variants are suspect (BE should detect them). PE-only hits in non-BE-coverable variants (e.g., transversions) are genuinely PE-unique.
