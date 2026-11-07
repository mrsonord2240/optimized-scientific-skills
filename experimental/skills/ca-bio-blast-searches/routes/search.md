# Submit a remote search

```bash
python scripts/blast_rid.py run --query query.fa --program blastn --db refseq_select_rna --hitlist 500 --expect 1e-10 --out hits.xml
```

Select blastn, blastp, blastx, tblastn or tblastx according to query and target types; see `references/programs.md`. For megablast pass `--megablast` with blastn. The RID CLI does not expose every qblast parameter: use `references/dna-patterns.md` for discontiguous search and `references/protein-patterns.md` for organism restriction or short-peptide scoring. A small hit-list limit may affect which hits are retained; retrieve a sufficiently large set and filter afterwards. Avoid concurrent remote jobs and do not poll more frequently than the source default 60 seconds.

Done when XML is saved and search program, parameters, database, date and query are recorded.
