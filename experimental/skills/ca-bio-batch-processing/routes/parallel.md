# Count FASTA files in parallel
```bash
python scripts/batch_cli.py parallel data '*.fasta' --workers 4
```
**Done when:** JSON includes a completed result for each discovered file in stable order. The worker operates on FASTA; it reads whole files in separate processes and streams records. The guarded main entry supports spawn platforms. See [parallel recipe](../references/operations.md#parallel-processing).
