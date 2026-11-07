# Search a local database

```bash
blastn -query query.fasta -db custom_db -task dc-megablast -word_size 11 -evalue 1e-10 -num_threads 8 -max_target_seqs 500 -outfmt '6 qseqid sseqid pident length qcovs qcovhsp evalue bitscore stitle' -out hits.tsv
```

Select task for the data: source recommends dc-megablast for cross-species DNA, megablast for high identity, and blastn-short for short queries; use blastp and a protein database for proteins. A limited target set can influence the retained hits, so use an appropriate sufficiently large limit before ranking. `qcovs` is subject-level query coverage; `qcovhsp` is per HSP. See `references/search.md` for masking, output and thread options.

Done when the result table is checked and task, database, thresholds and coverage measures are reported.
