# Multi-task modeling

```bash
chemprop train --data-path multitask.csv --task-type classification --target-columns CYP1A2_inhibition CYP2D6_inhibition CYP3A4_inhibition --save-dir multitask_model
```

Replace target columns with the dataset's actual endpoint names and verify current CLI flags. Compare single-task and multitask results on identical splits; see `references/multitask.md`.

Done when model and matched held-out comparisons are recorded.
