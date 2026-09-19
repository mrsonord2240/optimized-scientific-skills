# Single-Cell Data I/O - Usage Guide

## Overview

This skill covers reading, writing, creating, and converting single-cell objects across AnnData (Python), Seurat (R), and SingleCellExperiment (R). It emphasizes the decisions that prevent silent data loss: keeping the raw 10X matrix, choosing stable gene identifiers, picking a storage format, and converting between Python and R without dropping layers or transposing the matrix incorrectly.

## Prerequisites

See SKILL.md's Installation section for install commands (Python: scanpy/anndata/muon; R: Seurat, anndataR, zellkonverter, schard).

## Quick Start

Tell the AI agent what is needed:
- "Load the raw 10X matrix and keep the antibody-capture features"
- "Create an AnnData object from this count matrix and store raw counts"
- "Convert this h5ad to a Seurat object without losing the UMAP and layers"

## Example Prompts

### Loading Data
> "Read the raw_feature_bc_matrix folder using Ensembl gene IDs"

> "Load this Cell Ranger h5 and report cells x genes"

> "Load the 10X output but keep CRISPR guide and antibody features"

### Creating Objects
> "Build an AnnData from this matrix, put integer counts in a counts layer"

> "Create a Seurat v5 object from this sparse matrix with min.cells 3"

### Converting
> "Convert this AnnData to a SingleCellExperiment, keeping reducedDims and raw"

> "Move this Seurat object to h5ad for Python and verify no layers were dropped"

> "Why did my gene and cell axes swap after conversion?"

What the agent does and the traps it avoids (transpose, sparsity, lossy conversion, raw-vs-filtered) are in SKILL.md's Governing Principle, Common Errors, and API Defaults sections.

## Related Skills

- single-cell/preprocessing - QC, normalization, and HVG selection after loading
- single-cell/doublet-detection - per-sample doublet calling on raw counts after loading
- single-cell/clustering - dimensionality reduction and clustering on the loaded object
- single-cell/multimodal-integration - MuData/h5mu handling for CITE-seq and Multiome
- spatial-transcriptomics/spatial-data-io - SpatialData/zarr I/O for spatial omics
- workflows/scrnaseq-pipeline - end-to-end scRNA-seq pipeline that starts from data loading
