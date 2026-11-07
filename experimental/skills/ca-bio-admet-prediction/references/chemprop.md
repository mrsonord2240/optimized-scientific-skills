## chemprop D-MPNN for Custom Endpoints

When in-house data is available, train a target-specific model. chemprop provides a widely used open-source D-MPNN architecture with atom/bond features and optional molecular descriptors; benchmark it against appropriate baselines on the project's data.

**Goal:** Train a target-specific ADMET classifier or regressor on in-house bioassay data.

**Approach:** Use the installed Chemprop 2.x CLI with a scaffold split, a release-supported descriptor featurizer, replicated models, and an explicit prediction-time uncertainty/calibration workflow.

```python
# Chemprop 2.2 CLI; verify flags against the installed release.
# chemprop train --data-path data.csv --task-type classification \
#                --save-dir model_dir --split-type scaffold_balanced \
#                --molecule-featurizers v1_rdkit_2d_normalized \
#                --num-replicates 5 --ensemble-size 5
# chemprop predict --test-path test.csv --model-paths model_dir \
#                  --uncertainty-method ensemble --preds-path predictions.csv

# Or chemprop 2.x programmatic API (full programmatic API documented at chemprop.readthedocs.io)
# See chemoinformatics/qsar-modeling for the full chemprop 2.x training pipeline.
```

**Key:** Replicates and an ensemble estimator produce an uncertainty estimate, not automatic calibration. Fit and evaluate a documented calibrator on a separate calibration set when calibrated probabilities or intervals are required. Descriptor benefit must be demonstrated on the intended endpoint.
