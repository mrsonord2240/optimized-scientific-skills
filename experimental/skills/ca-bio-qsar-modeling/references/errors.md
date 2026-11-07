## Per-Tool Failure Modes

### Random split for QSAR

**Trigger:** Default sklearn `train_test_split`.

**Mechanism:** Compounds from same scaffold scatter across train/test; performance optimistic.

**Symptom:** Performance can drop from random splits to scaffold, time, external-series, or prospective evaluation; quantify the size rather than assuming it will be large.

**Fix:** Use `--split scaffold_balanced` in chemprop 2.x (or `--split_type scaffold_balanced` in chemprop 1.x legacy); or `scaffold_split` from `chemoinformatics/scaffold-analysis`.

### Class imbalance not handled

**Trigger:** 10:1 negative:positive ratio in dataset.

**Mechanism:** Default loss treats classes equally; model learns majority class.

**Symptom:** High accuracy but precision/recall on minority class poor.

**Fix:** Class-weighted loss; SMOTE; or report AUC/F1 not accuracy.

### Over-engineered features

**Trigger:** Including hundreds of descriptors (e.g., `rdkit_2d` not normalized).

**Mechanism:** Some descriptors dominate scaling; model overfits.

**Symptom:** Validation performance differs widely across runs; high feature importance noise.

**Fix:** In current chemprop 2.x use `rdkit_2d`, which is scaled by default, or supply a documented descriptor set with preprocessing fit only on the training data.

### Missing AD assessment

**Trigger:** Predicting on novel chemotypes without AD check.

**Mechanism:** Model extrapolates; predictions unreliable.

**Symptom:** On low-similarity compounds, explained variance or rank ordering can collapse even when MAE changes little because the novel subset has a narrower label range.

**Fix:** Predefine and validate one or more domain/uncertainty diagnostics, such as neighborhood similarity, ensemble disagreement, or conformal output. Stratify held-out R2 or Spearman correlation by the AD diagnostic (not MAE alone), and report what each diagnostic does and does not guarantee.

### chemprop 1.x vs 2.x confusion

**Trigger:** Code/tutorial from before late 2024.

**Mechanism:** Major API change: `chemprop_train` -> `chemprop train`; Python API redesigned.

**Symptom:** ImportError or different keyword arguments.

**Fix:** Use `chemprop --version`; check 2.x documentation; migrate APIs.

### Pretrained Transformer overhead without data benefit

**Trigger:** Adding a pretrained transformer without a matched baseline and deployment-relevant validation.

**Mechanism:** The pretrained representation, fine-tuning design, and endpoint may not provide additional transferable signal.

**Symptom:** No improvement over chemprop; slower training.

**Fix:** Compare against fingerprint and chemprop baselines on the same split, and retain the transformer only when the measured benefit justifies its cost.

### Validation leakage via standardization

**Trigger:** Standardization rules or learned preprocessing parameters are chosen or fit using validation/test data.

**Mechanism:** Test-set information influences representations, feature selection, scaling, or deduplication decisions.

**Symptom:** Re-fitting preprocessing on training data alone reduces held-out performance or changes membership across splits.

**Fix:** Freeze chemistry rules before evaluation and fit learned preprocessing on training data only. Apply the frozen pipeline to validation, test, and prospective compounds while preserving endpoint-relevant stereochemistry.

## Reconciliation: Classical RF vs chemprop vs Transformer

| Aspect | RF + ECFP4 | chemprop D-MPNN | MolFormer |
|--------|-----------|-----------------|-----------|
| Data regime | Useful baseline across sizes; especially important in small data | Compare when graph learning is plausible | Compare when pretrained representations and compute are justified |
| Interpretability | Fingerprint importance or SHAP, with bit-to-atom mapping caveats | Graph attribution requires a custom, validated implementation | Model-specific attribution requires validation |
| Uncertainty | Bootstrap or conformal wrapper | Ensemble disagreement; calibrate separately when needed | Method-dependent; validate empirically |
| Hardware | CPU | CPU or GPU depending on scale | Usually GPU for fine-tuning |
| OOD performance | Benchmark on the intended split/domain | Benchmark on the intended split/domain | Benchmark on the intended split/domain |
| Production deployment | Version-pinned sklearn artifact or service | Version-pinned native checkpoint/service; do not assume ONNX support | Version-pinned framework artifact/service |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| chemprop hangs at start | GPU OOM | Reduce batch_size; check CUDA |
| All predictions same value | Constant target | Standardize labels |
| AUC mismatched across folds | Split/data-order or PyTorch seed not fixed | `--data-seed 42 --pytorch-seed 42` (they control different randomness sources) |
| Test AUC = train AUC | No held-out data | Use scaffold_balanced split |
| Ensemble variance always small | Ensemble members insufficiently diverse | Check the documented seed behavior and training randomness for each replicate/member |
| SHAP fails on D-MPNN | Graph inputs are not compatible with the tree-model interface | Use a tested graph-attribution implementation or report the fingerprint baseline attribution |
| MolFormer fine-tune slow | All parameters trained | Use LoRA or freeze early layers |
| Calibration degrades held-out results | Calibrator overfit or distribution shifted | Refit on a proper calibration split and report uncalibrated and calibrated metrics |
