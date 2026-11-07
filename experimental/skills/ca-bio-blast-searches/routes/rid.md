# Resume a remote job

```bash
python scripts/blast_rid.py status RID
```

Use `python scripts/blast_rid.py submit --query query.fa --program blastn --db refseq_select_rna` to create an RID and `python scripts/blast_rid.py fetch RID --out hits.xml` after READY. Keep the RID; long WAITING times do not alone imply failure. The fetch writes the specified output path, so choose a fresh file. See `references/results-and-rid.md`.

Done when the RID is resolved or its failure or expiration is documented and XML saved if ready.
