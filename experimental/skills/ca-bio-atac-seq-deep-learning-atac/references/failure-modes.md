# Failure Modes and Troubleshooting

## chromBPNet -- Bias model mismatch

**Trigger:** Bias model trained or borrowed from a different sample (e.g. a K562 bias model on primary T cells).

**Mechanism:** The bias model captures Tn5 sequence preference, mostly cell-type-invariant, but chromatin context at cuts varies. Cross-sample bias models work with degraded fidelity.

**Symptom:** Footprints look right at CTCF but the nobias model still responds to Tn5 motifs, or TF-MoDISco of the nobias model returns Tn5-like or many GC-rich patterns.

**Fix:** Train the bias model on non-peak regions of the same library (`bias train`, raise `-b` if the bias model correlates negatively in peaks). If that is impossible, use an ENCODE `model.bias_scaled.*.h5` and accept the degradation; either way apply the QC gate.

## chromBPNet -- Insufficient training data

**Trigger:** Few deduplicated reads or peaks (uncited rule of thumb: < 50M reads, < 30k peaks).

**Mechanism:** Too few peaks give unstable gradients and few backgrounds for the joint loss; profile metrics such as JSD are read-depth sensitive.

**Fix:** Pool replicates; reduce capacity (`-fil`, `-dil`); or use the closest pre-trained ENCODE model and skip retraining.

## BPNet / chromBPNet -- DeepLIFT vs Integrated Gradients

**Trigger:** Computing per-base contributions for motif discovery.

**Mechanism:** chromBPNet uses DeepSHAP (DeepLIFT rescale rule with dinucleotide-shuffled references), which satisfies additivity; Integrated Gradients is stochastic and differs.

**Fix:** Use DeepLIFT/SHAP for TF-MoDISco; use IG only when DeepLIFT misbehaves on saturating activations. Record the choice and the head (counts vs profile).

## scBasset -- Cell embedding instability

**Trigger:** Few cells, very sparse cells, or few peaks.

**Mechanism:** Each cell has its own learned embedding row; sparse cells give noisy embeddings, and the preprocessing split needs >= 1,000 validation and test peaks.

**Fix:** Filter low-count cells (`min_counts`), keep enough peaks, and compare embeddings across seeds. For cluster-level questions train chromBPNet on pseudobulks instead.

## Enformer -- Reference tracks lack the target cell type

**Trigger:** Variant effects in a cell type absent from Enformer's training tracks.

**Mechanism:** Enformer predicts fixed reference tracks; a proxy track is only an approximation.

**Fix:** Choose a proxy track (`targets_human.txt`) and document it, or fine-tune (expensive).

## tangermeme -- Marginal effect vs in silico mutagenesis

**Trigger:** Asking for a "variant effect score" without a formula.

**Mechanism:** Marginal effect = ref vs alt at the SNP only. ISM = every base in the window mutated. Different magnitudes and questions.

**Fix:** State which. GWAS variants: marginal at the SNP. "Which bases matter": ISM.

## Common Errors and Solutions

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `FileNotFoundError` in `bias`/`train` for the `-n` file | No non-peak file; nothing else generates it | `chrombpnet prep nonpeaks -g -c -p -fl -o <prefix>` -> `<prefix>_negatives.bed` |
| `FileExistsError: <prefix>_auxiliary` | `prep nonpeaks` refuses to overwrite | Delete the directory or change the prefix |
| `bias pipeline` / `pipeline` never exits on CPU; log ends at `Generating profile shap scores` | Post-training DeepSHAP interpretation over all peaks | Use `bias train` / `train` (the shipped script) and run `contribs_bw` on a peak subset; models and metrics are written before that stage |
| `chrombpnet train` aborts with an assertion on the bias model's peak correlation | Bias model captured AT/GC content (pearsonr < -0.5) | Increase `-b` and retrain the bias model |
| `pred_bw` / `contribs_bw` `ValueError: cannot convert float NaN to integer` | Regions are not 10-column narrowPeak | Add the summit column (see method reference); also drop duplicate rows (`RuntimeError: entries ... out of order`) |
| variant-scorer `File has 4 columns but chrombpnet schema expects 5 columns` | 4-column variant list | `chr pos ref alt variant_id`, no header |
| variant-scorer `ModuleNotFoundError: pybedtools` | Undeclared dependency; pip build fails | `conda install -c bioconda pybedtools` |
| Variant effects ~5x too small | `log2` of a log-count output, or tuple output passed to tangermeme | `(y_alt - y_ref) / ln 2` on `CountWrapper(model)` output |
| tangermeme model has no `.sum` / returns lists | Raw `(profile, counts)` model output | Wrap with `CountWrapper` |
| `modisco motifs`: `Window (500) cannot be longer than the sequences` | Attributions are `(N, L, 4)` | Transpose to `(N, 4, L)` |
| `modisco report`: `tomtom executable could not be called` | MEME suite not on PATH | Put `tomtom` on PATH or use `-l` |
| modisco finds no motifs | Too few seqlets, or wrong contribution layout | Check shape `(N, 4, L)`; raise `-n`; use more peaks |
| Out of memory during training | Batch too large | `-bs` smaller (chromBPNet) |
| Predicted profile is constant | Model collapsed, or too few epochs | Train longer; check the peaks file is non-empty |
| Model does not converge | Peaks include chrM or blacklist regions | Pre-filter; chromBPNet does not filter them |
| Variant effects near zero | SNP outside the model's effective window | Centre the 2,114 bp window on the variant |
| Enformer wrong shape | Input not 196,608 bp | Pad/trim to exactly 196,608 bp |
| scBasset training: `ModelCheckpoint filepath must end in .weights.h5` | Keras 3 environment | Use the TF 2.15 / Keras 2 `scbasset` env |
| scBasset preprocess: `number sections must be larger than 0` | < 1,000 validation/test peaks | More peaks (> ~20k) or a lower cell-fraction filter |
| `python variant_scoring.py --help` fails | variant-scorer not cloned | Clone kundajelab/variant-scorer; run its `src/variant_scoring.py` with the `chrombpnet` env |
| `tangermeme` signature mismatch | API changed | `help(tangermeme.variant_effect.substitution_effect)` |
