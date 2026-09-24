# PTM Analysis - Usage Guide

## Overview
Localize and quantify post-translational modifications (phosphorylation, acetylation, ubiquitination, glycosylation) as three stacked inferences: which sites the enrichment chemistry actually captured, where on the peptide the modification sits (with its own false localization rate), and whether an abundance change survives subtracting the protein-level change. The central decision is never trusting a phospho-only fold-change as "regulation" without protein-level adjustment.

## Quick Start
Tell your AI agent what you want to do:
- "Load my MaxQuant Phospho (STY)Sites.txt, expand multiplicity, and keep class I sites"
- "Adjust phosphosite changes for protein abundance using MSstatsPTM and a paired global proteome"
- "Protein-adjust my TMT phosphoproteomics -- enriched and global runs are labelled plexes with a pooled reference channel"
- "Build a kinase-motif logo using an experiment-matched background, not the whole proteome"
- "Infer which kinases are active with KSEA from my site fold-changes"
- "Check whether my diGly sites are confounded by NEDD8/ISG15 or an iodoacetamide artifact"

## Example Prompts

### Site Identification and Localization
> "Filter Phospho (STY)Sites.txt to class I (localization probability >= 0.75) and report the residue distribution"

> "Expand the MaxQuant site-table multiplicity into Intensity___1/___2/___3 before any quantification"

> "Estimate an empirical global false localization rate for my phosphosites"

### Protein-Adjusted Quantification
> "Use MSstatsPTM groupComparisonPTM and call only ADJUSTED.Model hits as regulated"

> "Show me which apparent site changes are actually driven by protein abundance"

> "Compare phosphosite changes after drug treatment, adjusting for protein-level stabilization"

### Other PTMs
> "Treat my K-GG data as ubiquitin plus NEDD8 plus ISG15 and flag the chemistry confounds"

> "Map acetylation sites allowing four or more missed cleavages because acetyl-K blocks trypsin"

> "Disambiguate glycosite N->D from spontaneous deamidation"

### Motif and Kinase Activity
> "Run motif analysis with an experiment-matched S/T/Y background and render a sequence logo"

> "Run KSEA or PTM-SEA to infer active kinases and report z-scores with substrate counts"

## Where the detail lives
Decision tree, code, failure modes and install notes are all in `SKILL.md`; this guide only helps you choose the Skill and phrase the request.

## Related Skills
- peptide-identification - Identify modified peptides and run open/variable-mod search
- quantification - Underlying protein-level quant feeding the MSstatsPTM PROTEIN dataset
- differential-abundance - Moderated testing on the protein-level intensity matrix
- pathway-analysis/gsea - Enrichment scoring of regulated-site protein lists and PTM-SEA-style signatures
- data-visualization/sequence-logos - Render motif logos from the foreground/background windows
