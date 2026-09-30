---
name: bio-atac-seq-enhancer-gene-linking
description: Predict enhancer-gene regulatory connections from ATAC-seq using ABC, ENCODE-rE2G, HiChIP, or Cicero. Use when linking distal enhancers to target genes, choosing between contact-aware (ABC, ENCODE-rE2G), accessibility-only (ABC powerlaw, Cicero), and orthogonal (HiChIP H3K27ac, EpiMap) approaches, validating predictions against CRISPRi-FlowFISH gold-standard, or building cell-type-specific regulatory maps for fine-mapping or therapeutic target discovery.
license: MIT
category: Data Analysis
author: GPTomics
---

# Enhancer-gene linking

Predict which gene a distal accessible region regulates by combining accessibility activity, 3D contact frequency, and optionally sequence features. The result is a per-(enhancer, gene) score that can be thresholded. ABC and ENCODE-rE2G are the canonical predictors; ABC also runs without Hi-C by approximating contact with a distance powerlaw. Cicero reports co-accessibility from single-cell data, not contact or target genes. CRISPRi-FlowFISH is the experimental gold standard.

## Workflow

1. Inventory the inputs: ATAC (or DNase) and H3K27ac alignments, Hi-C/Micro-C (`.hic`, or a directory for `avg`/`juicebox`/`bedpe`) and its cell type, gene annotation, blacklist, genome build.
2. Choose the method from the table below. Match Hi-C to the cell type; if a proxy or fallback contact is used, name it and report the degradation.
3. For ABC, run [`scripts/run_abc.sh`](scripts/run_abc.sh): MACS2 summit peaks, candidate regions (summit +/- 250 bp, blocklist removed, TSS regions added so promoters are classified rather than dropped), quantile-normalised neighborhoods, prediction, and calibrated thresholding that removes non-self promoters. Set `ACCESS_TYPE` (`DHS` or `ATAC`) to match the accessibility input.
4. For ENCODE-rE2G, configure the Snakemake pipeline through `config/config.yaml` (`ABC_BIOSAMPLES` points to a biosamples TSV with accessibility, H3K27ac, `HiC_file`, `HiC_type`, `HiC_resolution`) and run `snakemake -j1 --use-conda`; there is no `cell_type=`/`atac_bw=` `--config` override. The model is auto-selected from the biosample inputs (see the method reference); do not pick one by tissue.
5. Threshold with the method's calibrated cut-off, cross-check methods with [`scripts/combine_predictions.py`](scripts/combine_predictions.py), and label single-method calls as exploratory.
6. Record method, tool versions, Hi-C source or fallback, thresholds, genome build, and validation.

## Method by available data

| Available data | Method |
|---|---|
| Accessibility + H3K27ac + matched Hi-C/Micro-C | ABC or ENCODE-rE2G (`HIC_TYPE=hic`) |
| Accessibility + H3K27ac, no cell-type Hi-C | ABC powerlaw contact (`HIC_TYPE=none`, the ABC-documented fallback), ABC average Hi-C (`avg`, needs the 58 GB ENCFF134PUN file), or ENCODE-rE2G with the megamap Hi-C; ENCODE-rE2G has no powerlaw or no-Hi-C model |
| Bulk accessibility only (no H3K27ac) | ABC with no `H3K27AC_BAM`, or ENCODE-rE2G `*_intact_hic` / `*_megamap` accessibility-only models; the ABC threshold table has separate rows |
| Single-cell ATAC | Pseudobulk by cluster (ABC recommends >= 3 million unique fragments, preferably 6 million) and run ABC; Cicero (atac-seq/co-accessibility) for co-accessibility hypotheses only |
| ATAC + H3K27ac HiChIP | FitHiChIP loops intersected with ABC |
| Multiome (ATAC + RNA, same cell) | LinkPeaks (Signac); SCENIC+ for TF networks |
| ENCODE cell type | Check for released ENCODE-rE2G predictions before running |
| Multi-cell-type scATAC | scBasset (atac-seq/deep-learning-atac) |
| Experimental validation wanted | CRISPRi-FlowFISH design, with predictions as hypotheses |

EpiMap, GeneHancer, and FANTOM5 are cell-type-agnostic baselines and sanity checks, never the primary call.

## Thresholds

Use the pipeline's calibrated cut-off; a single fixed value ignores input-dependent score distributions.

| Method | Threshold |
|---|---|
| ABC.Score / powerlaw.Score | Looked up by accessibility type, H3K27ac presence, and contact type in ABC's `reference/abc_thresholds.tsv` (0.012-0.027); `run_abc.sh` does this. 0.02 is the pipeline fallback when no row matches |
| ENCODE-rE2G.Score | The `threshold_*` file in the selected model directory (0.179-0.298; extended model 0.336), chosen for 70% recall on CRISPR-validated links |
| FitHiChIP loops | FDR < 0.05 and contact count >= 5 (convention; FDR < 0.01 and count >= 10 stringent) |
| Cicero co-accessibility | > 0.25 (> 0.5 stringent), convention |

For therapeutic target nomination require (a) a link passing its method's calibrated threshold, (b) agreement between two methods (ABC + ENCODE-rE2G, or ABC + HiChIP), and (c) experimental validation, CRISPRi-FlowFISH preferred. Report method intersections as high-confidence and unions as exploratory.

## Routed material

- [`scripts/run_abc.sh`](scripts/run_abc.sh): ABC pipeline configured by environment variables documented in its header.
- [`scripts/combine_predictions.py`](scripts/combine_predictions.py): intersects thresholded ABC and ENCODE-rE2G tables and flags HiChIP loop support.
- [`references/method-reference.md`](references/method-reference.md): method taxonomy, ABC mathematics, ENCODE-rE2G model selection, per-tool failure modes, reconciling disagreements, CRISPRi-FlowFISH validation, verified versions, common errors, and citations.
- [`references/usage-guide.md`](references/usage-guide.md): prerequisites, install commands, and request examples.

## Related skills

atac-seq/co-accessibility, atac-seq/atac-peak-calling, atac-seq/consensus-peakset, atac-seq/deep-learning-atac, atac-seq/single-cell-atac, hi-c-analysis/loop-calling, hi-c-analysis/contact-pairs, chip-seq/peak-calling, gene-regulatory-networks/scenic-regulons.

## License and provenance

This derived Skill retains GPTomics/bioSkills material by Domen Jemec, from commit `d91ed3d563019e649dc854c56ccd62551359488a`. See [`LICENSE`](LICENSE) for the preserved MIT notice.
