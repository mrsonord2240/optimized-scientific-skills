# PRIDICT2 batch interface
The source used a 2026-09 snapshot of PRIDICT2. Batch CSV columns are sequence_name and editseq; the latter contains edit notation. The default input directory is ./input. A fresh output directory must exist and contain no CSV files when --summarize is used. The source uses three cores and requires a model argument to --summarize. Inspect per-sequence output and logs: malformed input, no rows, or no successful designs can produce an empty summary despite exit zero. Actual output columns and score scale must match the chosen installed version.

The source reports a forced 5-prime G in Spacer-Sequence. This can differ from the corresponding genomic base; inspect the actual sequence and model version.
