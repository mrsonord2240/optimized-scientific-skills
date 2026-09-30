# Tool Decision Tree

Choose the tool and route by analytical goal. Commands, contracts and thresholds are in [`method-reference.md`](method-reference.md).

## By Goal

| Goal | Recommended approach |
|------|---------------------|
| Score ~100 GWAS SNPs for chromatin effects | Pre-trained ENCODE chromBPNet on the closest cell type; variant-scorer (null-based p-values) or the PyTorch route |
| Score one lead SNP at high resolution | PyTorch route (`substitution_effect`); a saturation-mutagenesis map with tangermeme was not exercised here |
| Identify TF binding motifs from a new cell type's ATAC | Train chromBPNet (QC gate), DeepLIFT/SHAP attributions, `modisco motifs`/`report` |
| Predict accessibility in a cell type not in training | Enformer (needs a matching reference track); chromBPNet trained on that cell type when data exist |
| Bias-correct a low-input ATAC library | chromBPNet corrected model (bias model from the same library) |
| Cell-type-specific enhancer prediction | chromBPNet per cell type; score candidate loci with the PyTorch route |
| Replace TOBIAS bias correction | Untested here; see SKILL.md step 8 |
| Per-cell scATAC embedding | scBasset on the cells x peaks matrix (`scbasset` env) |

## Quick Reference

**Variant effects:** use a pre-trained model if one exists, else train (`scripts/chrombpnet_pipeline.sh`); pass the QC gate; score with variant-scorer or `scripts/score_variants_torch.py`; for distal regulation beyond the model's 2,114 bp window add Enformer.

**Motif discovery:** train or load a model; QC gate; attributions -> `modisco motifs` -> `modisco report`; validate against JASPAR/HOCOMOCO.

**Bias correction:** `bias train` on non-peak regions; `train` with `-b bias.h5`; corrected signal from `pred_bw` for chosen regions.

## Version and API Stability Notes

- **chromBPNet:** `chrombpnet snp_score` is unavailable (use kundajelab/variant-scorer); `--help` lists the current subcommands.
- **tangermeme:** verify `substitution_effect` and `marginalize` with `help()`.
- **modisco-lite:** CLI only (`modisco motifs`, `modisco report`); replaces `kundajelab/tfmodisco`.
- **scBasset:** research code, not on PyPI; TF <= 2.15 for training.
- **Enformer:** `enformer-pytorch` (not in Kipoi); fine-tuning needs large GPUs and was not tested.
- **TensorFlow / PyTorch:** separate environments per tool; see SKILL.md Version Compatibility.
