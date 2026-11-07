# Attributions and motifs

```bash
python scripts/attributions_to_modisco_npz.py chrombpnet_nobias.h5 genome.fa peaks.narrowPeak motif_input 1500 counts
```

Apply QC first. Use the resulting one-hot and attribution arrays with modisco motifs and report as in `references/method-reference.md`. Arrays are N x 4 x 1000; choose counts or profile head deliberately and validate motif matches against a motif database.

Done when motif report exists and model and motif checks are reviewed.
