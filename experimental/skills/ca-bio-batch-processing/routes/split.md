# Split by record count
```bash
python scripts/batch_cli.py split input.fasta --format fasta --chunk-size 1000 --output split/chunk
```
**Done when:** reported chunks contain the expected record total and each has at most 1000 records. Create the parent directory first; choose a fresh output prefix because existing chunks can be overwritten. Chunk size must be a positive integer.
For safe prefix grouping use `python scripts/batch_cli.py prefix input.fasta --format fasta --output split-by-prefix`. It groups IDs before the first underscore, validates portable names, limits open handles to 64 and creates files exclusively. Inspect `split-prefix-manifest.json` for complete status, counts and errors; partial runs require review and a new directory. See [split recipes](../references/operations.md#split-files).
