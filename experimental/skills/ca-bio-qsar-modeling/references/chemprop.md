## chemprop 2.x Training (CLI)

**Goal:** Train five replicated chemprop runs, each containing a five-model D-MPNN ensemble with RDKit 2D descriptor features and a scaffold-balanced train/validation/test split.

**Approach:** For current chemprop 2.x, invoke `chemprop train` with `--molecule-featurizers rdkit_2d`, `--num-replicates 5`, `--ensemble-size 5`, and `--split scaffold_balanced`. Replicates repeat splitting/training with incremented seeds; they are not five-fold cross-validation. Confirm the exact flags with `chemprop train --help` because the v2 CLI continues to evolve.

```bash
# chemprop 2.x CLI (current): use 'chemprop train' (space; dashes not underscores)
chemprop train \
    --data-path data.csv \
    --task-type classification \
    --save-dir model_dir \
    --molecule-featurizers rdkit_2d \
    --num-replicates 5 \
    --ensemble-size 5 \
    --epochs 50 \
    --batch-size 128 \
    --split scaffold_balanced \
    --split-sizes 0.8 0.1 0.1 \
    --metric roc \
    --data-seed 42 \
    --pytorch-seed 42

# chemprop 1.x legacy CLI (for backwards reference):
# chemprop_train --data_path data.csv --dataset_type classification ...
```

Key flags (chemprop 2.x):
- `--molecule-featurizers rdkit_2d`: include current v2 RDKit descriptors, which are scaled by default (the legacy v1-normalized generator is `v1_rdkit_2d_normalized`)
- `--num-replicates 5`: repeat the split/training workflow with successive seeds; this replaced `--num-folds` in chemprop 2.1
- `--ensemble-size 5`: train five models per replicate for an ensemble prediction
- `--split scaffold_balanced`: prevent scaffold leakage (was `--split_type` in 1.x)
- `--split-sizes 0.8 0.1 0.1`: 80/10/10 train/val/test
- `--data-seed 42` and `--pytorch-seed 42`: respectively fix split/data-order randomness and PyTorch training randomness; set both when reproducing a run

Total models: 25 (5 replicates x 5 ensemble members). Report which predictions are aggregated. `chemprop predict --uncertainty-method ensemble` writes `pred_*_unc` as ensemble **variance**, not standard deviation; take its square root only when a standard-deviation diagnostic is required. Neither is a calibrated guarantee.

At prediction time, uncertainty output is opt-in and requires the actual saved model paths:

```bash
chemprop predict --test-path test.csv \
    --model-paths path/to/model_1.ckpt path/to/model_2.ckpt \
    --molecule-featurizers rdkit_2d \
    --uncertainty-method ensemble \
    --preds-path predictions.csv
```

Every molecule featurizer used at training must be supplied again at prediction. If it is omitted here, a model trained with `rdkit_2d` can fail with a 517-versus-300 tensor-shape error because prediction supplied graph features without the 217 descriptor features.
