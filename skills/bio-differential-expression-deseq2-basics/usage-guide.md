# DESeq2 Basics - Usage Guide

## Overview

This skill covers differential expression analysis using DESeq2, the most widely-used Bioconductor package for analyzing RNA-seq count data. DESeq2 uses a negative binomial model to test for differential expression between experimental conditions.

## Quick Start

Tell your AI agent what you want to do:
- "Run DESeq2 on my count matrix with treated vs control comparison"
- "Analyze differential expression controlling for batch effects"
- "Get significantly differentially expressed genes with padj < 0.05"

## Example Prompts

### Basic Analysis
> "Create a DESeqDataSet from my count matrix and sample metadata"

> "Run the standard DESeq2 workflow on this RNA-seq data"

> "Apply log fold change shrinkage to my DESeq2 results"

### Design Formulas
> "Set up DESeq2 with batch correction"

> "Create an interaction model for genotype and treatment"

> "Compare multiple conditions against a control"

### Results
> "Extract genes with adjusted p-value < 0.05 and |log2FC| > 1"

> "Order results by significance"

> "Export normalized counts and DE results"

## What the Agent Will Do

It follows the Standard Workflow in `SKILL.md`; install notes are under Version Compatibility there.

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Count matrix | Integer matrix | Genes (rows) x Samples (columns) |
| Sample metadata | Data frame | Rownames matching count column names |
| Design formula | R formula | Variables from metadata (~condition) |

The tips (reference level, explicit `name=`/`contrast=`, shrinkage, LRT, pseudobulk, prokaryotic caveats) live in `SKILL.md`.

## Related Skills

- edger-basics - Alternative DE analysis with edgeR
- de-visualization - MA plots, volcano plots, heatmaps
- de-results - Extract and export significant genes
- pathway-analysis/go-enrichment - Functional analysis of DE gene lists
- expression-matrix/counts-ingest - Loading count matrices from various formats
