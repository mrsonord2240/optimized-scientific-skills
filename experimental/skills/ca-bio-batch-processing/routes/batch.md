# Count a batch
```bash
python scripts/batch_cli.py count data '*.fasta' --format fasta
```
**Done when:** each discovered input is listed with a streamed count in stable order; check the pattern and expected inventory. Run from this skill directory. Biopython 1.83+ is the source target; core file recipes read plain text, not gzip.
Use `python scripts/batch_cli.py merge data '*.fasta' --format fasta --output merged.fasta` to stream a merge; decide how to handle duplicate IDs before indexing. `python scripts/batch_cli.py summary data '*.fasta' --format fasta --output summary.csv` writes count and length statistics, including a header-only file for no matches. `python scripts/batch_cli.py convert genbank '*.gb' --format genbank --output fasta --output-format fasta` creates one file per input; GenBank to FASTA drops features, annotations and qualifiers. Choose new paths: merge, conversion and summary may overwrite existing files. See [operation recipes](../references/operations.md), including source-tracking merge.
