# Batch shortcut check

```bash
python scripts/batch_checks.py data.csv --label label --batch-col batch --positive AML --out batch_results.tsv
```

Needs a batch column (two or more batches, at least 2 samples each). Writes a summary and a per-batch table: batch predictability (cross-validated one-vs-rest AUC for predicting the batch from the features), label-vs-batch chi-square p, random-split label AUC, leave-one-batch-out per-batch AUCs, their mean and the pooled AUC. It prints a red flag when the batch is predictable (AUC above 0.9) and associated with the label.

- Read the leave-one-batch-out per-batch mean, not the pooled AUC; pooling rewards a batch-level score shift that tracks each batch's case rate.
- A held-out batch with one class has no AUC (NaN, noted). A single leave-one-batch-out mean is noisy (spread about 0.1 at 2, 3 and 6 batches), so one run is not a verdict.
- If batch is confounded with the outcome, no correction rescues the design: ComBat cannot, and nested CV is still optimistic. Fix it at the design stage by balancing batches across outcome.
- Batch correction, if any, goes inside the training fold.

Done when `batch_results.tsv` exists and the answer states batch predictability, the label-batch association and the leave-one-batch-out per-batch mean.
