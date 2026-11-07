# Single-cell embeddings

```bash
python scBasset/bin/scbasset_preprocess.py --ad_file data.h5ad --input_fasta genome.fa --out_path processed
```

Use the compatible scBasset environment and preprocessing/training steps in `references/method-reference.md`. Input h5ad needs chr, start and end peak columns; verify adequate validation and test peaks and compare across seeds.

Done when embeddings and model quality are reviewed.
