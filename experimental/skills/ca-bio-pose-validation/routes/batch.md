# Batch pose QC

```bash
python scripts/pose_qc_batch.py receptor.pdb poses1.sdf poses2.sdf > batch_qc.txt
```

The batch helper reports the first PB-valid pose per file in input rank order. Inspect failures and expected checks separately; use `references/interpretation.md` for RMSD comparison and troubleshooting.

Done when retained poses and excluded or unscored inputs are reviewed.
