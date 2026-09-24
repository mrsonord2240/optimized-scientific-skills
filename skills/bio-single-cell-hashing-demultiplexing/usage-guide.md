# Hashtag Demultiplexing and Cross-Sample Doublet Calling - Usage Guide

## Overview

Hashtag demultiplexing assigns pooled, hashed cells back to their sample of origin and calls cross-sample doublets directly from the HTO count matrix. Each sample is labeled before pooling with a unique oligo-tagged reagent (CITE-seq antibody HTO, MULTI-seq lipid/cholesterol tag, or CellPlex CMO), so a singlet is dominated by one tag above background and a cross-sample doublet shows two tags high. This skill covers Seurat HTODemux and MULTIseqDemux, scanpy hashsolo, pegasus/demuxEM, GMM-Demux, and demuxmix, when to pick each, and how hashtag demux relates to genetic demux and expression-based doublet detection.

## Quick Start

Tell your AI agent what you want to do:
- "Assign my hashed cells back to their sample of origin"
- "Call cross-sample doublets from my HTO counts"
- "Run HTODemux on my Seurat HTO assay"
- "Demultiplex my MULTI-seq lipid-tagged samples"
- "I have lots of Negatives and weak staining - which demux method should I use?"

## Example Prompts

### Hashtag demultiplexing
> "Normalize my HTO assay with CLR and run HTODemux to assign samples and flag doublets"
> "Demultiplex these MULTI-seq tags with MULTIseqDemux using automated thresholding"
> "Use hashsolo in scanpy to assign samples when I only have two hashtags"

### Robustness to ambient and bad staining
> "My HTODemux call gives a huge Negative pile - rescue it with a method that models ambient background"
> "Run demuxEM on my nucleus-hashing data using empty droplets to estimate background"
> "Use demuxmix with the detected-gene count to handle uneven staining across tags"

### Choosing a modality
> "I did not hash but my samples are different donors - should I use genetic demultiplexing instead?"
> "Combine my hashtag doublet calls with expression-based doublet detection"

## What the Agent Will Do

See SKILL.md: Choosing a demultiplexing modality, Choosing a hashtag caller, the per-method sections, Common Errors and Related Skills. Install commands are in SKILL.md's Install section.
