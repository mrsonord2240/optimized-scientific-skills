# Prime-Editing Screens - Usage Guide

## Overview

Decision-grade design and analysis of pooled prime-editor screens. Covers pegRNA design with PRIDICT (Mathis 2023 Nat Biotechnol 41:1151) and PRIDICT2 (Mathis 2025 Nat Biotechnol 43:712, online June 2024); PE chemistry selection (PE2, PE3, PEmax, PE5max); pegRNA architecture (spacer + scaffold + PBS + RTT); chromatin context as a major locus-level determinant; PRIME pooled-screen methodology (Ren 2023); MOSAIC in situ saturation mutagenesis; CRISPResso2 quantification; and cross-modal validation with base-editor screens.

## Prerequisites

```bash
# PRIDICT2 (web + CLI)
git clone https://github.com/uzh-dqbm-cmi/PRIDICT2
# OR use web interface https://pridict.it/

# CRISPResso2 for amplicon analysis
conda install -c bioconda crispresso2   # not on PyPI
# or, as SKILL.md's block assumes, the container: pinellolab/crispresso2:latest

# ePRIDICT for the chromatin-context check (Linux/macOS only -- pybigwig has no Windows wheel;
# use WSL on Windows). Needs the model's ENCODE bigWigs too: ~5.3 GB for the light model.
git clone https://github.com/Schwank-Lab/epridict
conda create -n epridict -c conda-forge python=3.10 numpy pandas scipy tqdm joblib xgboost pybigwig
# then: ./epridict_download_encode.sh light   (see references/epridict-chromatin.md)

# Helpers (biopython is not used -- PBS/RTT extraction is a hand-rolled reverse complement)
pip install pandas numpy
```

Required inputs:
- Intended variant list with genomic coordinates
- Target genome reference (e.g. GRCh38)
- Cell line for PE chemistry validation
- The ePRIDICT model's ENCODE chromatin tracks for chromatin-aware filtering (downloaded with the
  tool, not a user-supplied ATAC-seq track)

## Quick Start

Tell the AI agent what to do:
- "Design a pegRNA library for 320 specific ClinVar MLH1, MSH2, MSH6, PMS2 variants. Use PRIDICT2 to score; select top 3 pegRNAs per variant"
- "Compare BE vs PE for installing a C->T variant at TP53 R273H. Recommend chemistry based on editing window + bystander"
- "Diagnose: my PE screen has 80% scaffold incorporation. Re-design pegRNAs with longer RTT and verify with PRIDICT2"
- "Run CRISPResso2 in prime-editor mode on amplicon FASTQ. Report intended-edit %, scaffold incorporation, indel byproduct per pegRNA"
- "Cross-validate PE screen hits with parallel BE screen at the same variants; identify PE-only hits in transversion variants"

## Example Prompts

### pegRNA Library Design

> "Build pegRNA library for 500 ClinVar MMR gene variants (MLH1, MSH2, MSH6, PMS2). Score each candidate pegRNA with PRIDICT2. Select 3 best pegRNAs per variant with predicted efficiency >50%. Output peg_library.csv with spacer, PBS, RTT, predicted_efficiency."

> "For installing MLH1 c.677A>G as a single intended variant: PE (no bystander) vs BE (no PAM-restriction). Run PRIDICT2 on PE candidates; check editing window for BE; recommend best chemistry."

### Screen Analysis

> "Run CRISPResso2 in prime-editor mode on my 50-pegRNA pilot. Report intended_edit_pct, scaffold_incorp_pct, indel_pct per pegRNA. Flag pegRNAs with intended_edit <5% for filtering or re-design."

> "Apply MAGeCK MLE on my PE screen counts after filtering to efficiency-passing pegRNAs. Aggregate per-pegRNA LFCs to per-variant scores."

### Chromatin-Aware Design

> "Score the loci my PRIDICT2-selected pegRNAs target with ePRIDICT (K562). Flag pegRNAs whose ePRIDICT percentile is low despite a high PRIDICT2 score for empirical pilot validation before library order. Mathis 2025 documented chromatin as a major locus-level determinant beyond PRIDICT sequence prediction, which is why ePRIDICT is meant to be combined with PRIDICT2.0."

### MOSAIC / Saturation

> "Design MOSAIC saturation library tiling TP53 DNA-binding domain (amino acids 94-294). For each residue, design pegRNAs covering all 19 possible amino acid changes (or at least all functionally relevant). Include variants from ClinVar."

### Cross-Validation

> "I have parallel BE and PE screens at the same 200 variants. Intersect hits at FDR <0.05; identify high-confidence variants (both methods, same direction). Flag PE-only hits in non-BE-coverable variants (transversions) as genuine PE-unique."

### Diagnostics

> "PE screen has 30% scaffold incorporation across library. Diagnose: RTT too short relative to PBS, or RT processivity issue with this cell line? Recommend pegRNA re-design or chemistry switch."

> "PRIDICT2 predicts 65% efficiency for my top pegRNAs, but pilot CRISPResso2 shows 8%. Score those loci with ePRIDICT to test chromatin context as the cause."

## What the Agent Will Do

1. Confirm intended variants are PE-installable (NGG PAM within 30 nt; achievable transversion/multi-base if BE-uncoverable)
2. Decide chemistry: PE2 standard; PEmax for high-efficiency; PE3 with second nick if expression supports
3. Design pegRNA candidates: spacer + PBS (11-13 nt, 40-55% GC) + RTT (10-20 nt) + scaffold
4. Run PRIDICT2 to predict per-pegRNA efficiency, indel rate, scaffold incorporation
5. Filter library to top 3 pegRNAs per variant with predicted efficiency >50%
6. Score the selected targets' chromatin context with ePRIDICT; flag low-percentile loci for pilot
7. Pilot the library at 20-50 representative loci; verify empirical efficiency matches PRIDICT2
8. Scale to full screen at MOI 0.3 in PE-expressing cell line
9. Endpoint amplicon sequencing of each target locus
10. CRISPResso2 quantification: intended-edit %, scaffold incorp %, indel %
11. Filter library to >5% intended-edit pegRNAs
12. Per-variant hit calling with MAGeCK MLE / drugZ
13. Cross-validate with BE screen if applicable
14. Output: per-variant fitness, BE-PE concordance, validation strategy

## Tips

- PRIDICT2 sequence prediction is a starting point, not a guarantee. Chromatin context is a major locus-level determinant beyond sequence (pair PRIDICT2.0 with ePRIDICT); always pilot at representative sites.
- For high-stakes variant-function calls, validate with both BE (where applicable) and PE. Concordant hits are high confidence; PE-only at BE-coverable variants are suspect.
- PE is the gold standard for transversions (C->G, C->A, G->C, G->A) where BE doesn't work; for these variants, PE is the only choice.
- Scaffold incorporation is a design issue, not a chemistry issue. Always check post-PE; if >5%, re-design pegRNAs with longer RTT.
- PE3 adds a nicked sgRNA on the opposite strand; this increases editing efficiency but doubles the indel risk. Use PE2 unless PEmax + PRIDICT2 doesn't reach target efficiency.
- For genome-scale saturation editing, MOSAIC (2024 preprint, demonstrated on BCR-ABL1 and the IRF1 UTR) is a high-throughput prime-editing option: 1000s of variants per cell line; smaller per-variant cell counts; suitable for drug-resistance variant scanning. For established saturation genome editing of BRCA1, the reference approach remains SGE (Findlay 2018).
- Pilot PRIDICT2 predictions in your specific cell line. The deep-learning model was trained on HEK293T (MMR-deficient) and K562 (MMR-proficient); other lines may have different chromatin and RT-activity profiles.
- For BE-PE comparison, design parallel libraries targeting the same variants. Concordant hits represent strong biology.

For the BE/PE/Cas9 decision tree, see SKILL.md's "Cas9 vs BE vs PE for Variant Installation".
For efficiency/QC thresholds, see SKILL.md's "Quantitative Thresholds".

## Validation Checklist

- [ ] Pre-synthesis PRIDICT2 filter (>50%)
- [ ] Chromatin context cross-referenced (ePRIDICT on the PRIDICT2-selected positions)
- [ ] Pilot empirical efficiency match >70% of PRIDICT2 prediction
- [ ] CRISPResso2 PE mode for quantification
- [ ] Scaffold incorporation <5% library-wide
- [ ] Per-variant 3+ pegRNAs in library
- [ ] BE cross-validation (if BE-coverable variants)

## Related Skills

- crispr-screens/library-design - pegRNA library design
- crispr-screens/base-editing-analysis - Orthogonal BE for variant attribution
- crispr-screens/crispresso-editing - CRISPResso2 PE mode
- crispr-screens/hit-calling - Per-variant hit calling
- crispr-screens/screen-qc - Editing-efficiency QC
- variant-calling/variant-annotation - Annotate edited variants
- clinical-databases/clinvar-lookup - Pathogenicity annotation
