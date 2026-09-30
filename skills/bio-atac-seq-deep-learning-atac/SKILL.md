---
name: bio-atac-seq-deep-learning-atac
description: Sequence-based deep learning for ATAC-seq using chromBPNet, BPNet, scBasset, or Enformer. Use when correcting Tn5 bias with neural networks beyond k-mer models, predicting per-base accessibility profiles, scoring in silico variant effects at GWAS or rare-variant SNPs, discovering motifs via DeepLIFT/TF-MoDISco from a trained model, or generating cell-type-specific accessibility predictions for unobserved cell states.
tool_type: python
primary_tool: chrombpnet
license: MIT
category: Data Analysis
author: GPTomics
---

## Version Compatibility

Verified 2026-09-30. No single environment satisfies every tool; keep three:

| Env | Packages | Used for |
|---|---|---|
| `chrombpnet` (Python 3.8) | chrombpnet 1.0.1 (pins tensorflow 2.8.0, keras 2.8.0, numpy 1.23.4), kundajelab/variant-scorer (git 0e1e341), pybedtools (conda; the pip build fails), bedtools and UCSC `bedGraphToBigWig` (conda `bedtools`, `ucsc-bedgraphtobigwig`), MEME suite (`tomtom`) | training, QC, `pred_bw`, `contribs_bw`, variant-scorer |
| `torch` (Python 3.11) | torch >= 2.1 (cu128 build for RTX 50-series), bpnet-lite 1.0.0, tangermeme 1.5.0, modisco-lite 2.4.0, enformer-pytorch 0.8.12, pyfaidx | attributions, TF-MoDISco, variant effect, Enformer |
| `scbasset` (Python 3.11) | tensorflow-cpu 2.15.1, keras 2.15.0, numpy 1.26.4, pandas < 3, anndata 0.11.4, psutil, scBasset (git aed3a6f, not on PyPI) | scBasset only |

TensorFlow 2.8 has no CUDA 12 build, so chromBPNet trains on CPU (the GPU helps only in the `torch` env). The PyPI name of TF-MoDISco-lite is `modisco-lite`. Deep-learning tooling evolves rapidly; verify with `pip show <package>` and `<tool> -h`.

## Workflow

1. **Choose the tool** by goal ([`references/tool-decision-tree.md`](references/tool-decision-tree.md)). Most common: chromBPNet for Tn5 bias correction and variant effects; Enformer for long-range distal regulation.
2. **Train** (`chrombpnet` env): [`scripts/chrombpnet_pipeline.sh`](scripts/chrombpnet_pipeline.sh) runs splits, non-peak generation, bias model, accessibility model and QC, and optionally variant scoring. It needs a 10-column narrowPeak and generates the non-peaks itself. Contracts and CPU runtimes: [`references/method-reference.md`](references/method-reference.md).
3. **Or use a pre-trained ENCODE model** (Keras `.h5`; e.g. GM12878 ATAC ENCSR095QNB, models tarball ENCFF935USI, `fold_*/model.chrombpnet_nobias.*.h5`). It loads in the `chrombpnet` env and in PyTorch through bpnet-lite. Match cell type and the 2,114 bp input window.
4. **Gate on QC** ([`scripts/chrombpnet_qc_gate.py`](scripts/chrombpnet_qc_gate.py)): do not score variants or interpret motifs from a model that fails it.
5. **Variant effects**: variant-scorer (5-column variant list, null-based p-values) or the PyTorch route ([`scripts/score_variants_torch.py`](scripts/score_variants_torch.py)). Effects are log2 count fold changes computed on the log-count head; see the method reference for both routes and the calling rule.
6. **Motifs**: DeepLIFT/SHAP attributions ([`scripts/attributions_to_modisco_npz.py`](scripts/attributions_to_modisco_npz.py) or `chrombpnet contribs_bw`), then `modisco motifs` and `modisco report`, then validate against JASPAR/HOCOMOCO.
7. **Enformer / scBasset**: [`scripts/enformer_variant_effect.py`](scripts/enformer_variant_effect.py); scBasset commands are in the method reference.
8. **Footprinting integration** (untested here): the chromBPNet corrected track can replace ATACorrect output for TOBIAS; `pred_bw` writes signal only for the regions you request.

### When to Use This Skill

**Yes, use deep learning for:**
- Variant effect prediction at GWAS causal SNPs (chromBPNet, or Enformer for distal)
- De novo motif discovery beyond MEME/HOMER (chromBPNet + DeepLIFT + TF-MoDISco)
- Tn5 bias correction at weak TF motifs or low-complexity sequence
- Per-cell accessibility embeddings in scATAC (scBasset)
- Cross-cell-type prediction for states without data (Enformer)

**No, use classical pipelines for:**
- Peak calling (MACS3), differential accessibility (DiffBind/DESeq2)
- Standard motif discovery in complex backgrounds (MEME/HOMER)
- Bulk TF activity (chromVAR at cluster level)

## Reference Routing

- [`references/failure-modes.md`](references/failure-modes.md): symptoms, causes and fixes, including the CLI errors listed there.
- [`references/method-reference.md`](references/method-reference.md): architecture comparison, pipeline contracts and runtimes, QC thresholds, variant-effect and TF-MoDISco routes, Enformer, scBasset, reconciliation of conflicting predictions.
- Full-scale chromBPNet training and interpretation, Borzoi, and scBasset on TF >= 2.16 were not run. Use pre-trained ENCODE models when the cell type exists; train only when it does not.

## Related Skills

- atac-seq/atac-peak-calling — Classical peak calling input
- atac-seq/footprinting — chromBPNet bias correction as a TOBIAS alternative
- atac-seq/motif-deviation — chromVAR vs scBasset for per-cell motif activity
- atac-seq/single-cell-atac — scBasset integration with sc workflow
- atac-seq/enhancer-gene-linking — Variant effect feeds enhancer scoring
- atac-seq/allele-specific-accessibility — DL-predicted effects vs observed allelic imbalance
- causal-genomics/fine-mapping — Downstream use of variant effect scores
- machine-learning/biomarker-discovery — General ML patterns
- gene-regulatory-networks/scenic-regulons — Motif discovery + TF networks
