# Long-Read Splicing - Usage Guide

## Overview
Alternative splicing from PacBio Iso-Seq (HiFi, Kinnex/MAS-Iso-seq) and Oxford Nanopore (direct cDNA, direct RNA) long reads: minimap2 alignment, isoform discovery and quantification (FLAIR, IsoQuant, Bambu), classification (SQANTI3), event-level testing (rMATS-long) and DTU. Commands, decision tree, install notes and failure modes are in `SKILL.md`. For a parameterized end-to-end starting point, see [`examples/longread_splicing_pipeline.sh`](examples/longread_splicing_pipeline.sh).

## Quick Start
Tell your AI agent what you want to do:
- "Discover and quantify isoforms from my PacBio Iso-Seq HiFi data"
- "Run FLAIR pipeline (correct -> collapse -> quantify -> diffSplice) on ONT direct cDNA"
- "Check whether my 10-nt microexon survives alignment, and rescue it if not"
- "Use Bambu for annotation-aware discovery + quantification with calibrated NDR"
- "Process MAS-Iso-seq + 10X 5' single-cell long-read data for full isoform per cell"

## Example Prompts

### PacBio Iso-Seq Discovery
> "Align HiFi reads with minimap2 -ax splice:hq, run IsoQuant for de-novo discovery, classify with SQANTI3."

### ONT Direct cDNA (unstranded by default; PCS-114, PCB-114)
> "Align ONT direct cDNA with minimap2 -ax splice -k14 (no -uf: the reads are unoriented), run FLAIR correct->collapse->quantify->diffSplice between conditions."

### Discovery + Quant in One
> "Use Bambu in R for joint isoform discovery and quantification with NDR=0.1 for balanced novel discovery."

### Event-Level on Long Reads
> "Use rMATS-long to identify SE/A5SS/A3SS/MXE/RI events from long-read isoform GTFs between two groups."

### Single-Cell Long-Read
> "Demultiplex MAS-Iso-seq Kinnex arrays with skera, then run lima -> isoseq refine for per-molecule FLNC reads."

### DTU on Long-Read Counts
> "Run DRIMSeq + stageR DTU on my FLAIR quantify counts matrix."

### Microexon-Specific
> "Confirm microexon inclusion in neural samples using long-read coverage; align with the annotation junctions and cross-validate with short-read junctions."

## Related Skills

- splicing-quantification - Short-read PSI for cross-validation
- isoform-switching - DTU framework on long-read counts
- single-cell-splicing - MAS-Iso-seq + 10X integration
- long-read-sequencing/isoseq-analysis - PacBio Iso-Seq general pipeline
- long-read-sequencing/long-read-alignment - minimap2 splice:hq details
- long-read-sequencing/long-read-qc - QC for long-read data
