# In-house endpoint models

```bash
chemprop train --data-path data.csv --task-type classification --save-dir model_dir --split-type scaffold_balanced --molecule-featurizers v1_rdkit_2d_normalized --num-replicates 5 --ensemble-size 5
```

Check installed 2.x flags before adapting the training and prediction recipe in `references/chemprop.md`. Use deployment-relevant splits, appropriate baselines, a separate calibration set, and documented applicability-domain assessment.

Done when model evaluation, uncertainty and calibration are reported for the intended data.
