# HMM or mapped columns

```bash
trimal -in input.fasta -out trimmed.fasta -gappyout -colnumbering > kept_columns.txt
```

`-gappyout` is the source choice for HMM preparation. The mapping gives zero-based original columns; ClipKIT logs use one-based positions. For alternative trimAl criteria and sequence-overlap filtering, read `references/trimal.md`. Record the version and explicit mode; `-automated1` can choose different modes across releases. Compare input and output sequence IDs if using overlap filters.

Done when the output and retained column mapping are checked and the selection mode recorded.
