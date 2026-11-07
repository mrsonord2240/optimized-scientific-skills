# Filter alignment

```bash
python scripts/filter_sequences.py input.fasta output.fasta --max-gap-fraction 0.1
```

Choose `--id-pattern`, `--max-gap-fraction`, and `--dedup` as needed; format can be set with `--format`. A filter removing every row fails. For column filtering, coordinate retention and reliability masking, use `references/trimming-and-reliability.md`; for row-filter details use `references/sequence-filtering.md`.

To remove columns with at least 50% gaps and then rows above 20% gaps, using source defaults:

```bash
python scripts/clean_columns.py input.fasta cleaned.fasta
```

Done when the requested filtered alignment is written and retained coordinates are reported.
