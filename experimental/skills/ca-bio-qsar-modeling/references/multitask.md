## Multi-Task QSAR

Train multiple related endpoints jointly:

```python
df = pd.DataFrame({
    'smiles': [...],
    'CYP1A2_inhibition': [...],
    'CYP2D6_inhibition': [...],
    'CYP3A4_inhibition': [...],
})
df.to_csv('multitask.csv', index=False)
```

```bash
chemprop train --data-path multitask.csv --task-type classification \
               --target-columns CYP1A2_inhibition CYP2D6_inhibition CYP3A4_inhibition \
               --save-dir multitask_model
```

Multitask learning can help when endpoints share predictive signal or data, but negative transfer is also possible. Compare single-task and multitask models under identical splits rather than assuming improvement from endpoint relatedness.
