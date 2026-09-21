# Isoform Switching - Usage Guide

## Overview
Identify shifts in which transcript a gene predominantly uses between conditions (differential transcript usage, DTU), and predict what the switch does to the protein (NMD sensitivity, ORF, domains, coding potential). This is not differential gene expression (DGE) or differential transcript expression (DTE). Tools, install commands, thresholds, the count route, and failure modes are in `SKILL.md`; `examples/isoform_switch_analysis.R` runs the whole workflow (built-in demo data, or your own Salmon output, GTF, FASTA and sample table).

## Quick Start
Tell your AI agent what you want to do:
- "Identify isoform switches between conditions and predict NMD/domain consequences"
- "Run the DTU pipeline with proper stageR multi-stage FDR control"
- "Use swish on Salmon Gibbs samples to incorporate quantification uncertainty"
- "Find genes where a poison exon switch reduces functional protein"

## Example Prompts

### IsoformSwitchAnalyzeR Workflow
> "Import my Salmon quantification (sample table attached), pre-filter for expressed isoforms, test for switches with a dIF cutoff of 0.1, and annotate consequences with whatever annotators I have (CPC2, Pfam, ...)."

> "Generate switch plots for the top 25 hits and compute global consequence enrichment."

### Manual DTU Pipeline
> "Run DRIMSeq filter -> DEXSeq exon-bin test -> stageR two-stage FDR for proper gene+transcript-level multiple testing."

### Inferential-Uncertainty-Aware
> "Re-run Salmon with --numGibbsSamples 20, then use fishpond/swish for DTE that propagates quantification uncertainty."

### NMD-Focused
> "Find isoform switches where the alternative form is NMD-sensitive (PTC >50nt upstream of last exon-exon junction)."

> "Identify poison-exon switches in SR/hnRNP genes that autoregulate via AS-NMD."

### Long-Read Integration
> "Run IsoformSwitchAnalyzeR on PacBio Iso-Seq transcript counts to bypass quantification uncertainty."

## Related Skills

- differential-splicing - Event-level (rMATS, leafcutter, MAJIQ) complement to DTU
- splicing-quantification - PSI is a 1D projection of DTU shifts
- splice-variant-prediction - Connects SpliceAI predictions to specific isoforms
- long-read-splicing - Full-isoform DTU bypasses transcript-quant uncertainty
- pathway-analysis/go-enrichment - Pathway enrichment of switching genes
- rna-quantification/alignment-free-quant - Salmon with `--numGibbsSamples` for swish
