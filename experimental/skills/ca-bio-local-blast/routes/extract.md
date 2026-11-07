# Extract or taxonomically restrict hits

```bash
blastdbcmd -db custom_db -entry_batch hit_accessions.txt -out hits.fasta
```

Supply one requested subject accession per line; indexed identifiers require `-parse_seqids` at build time. For `-taxids` or `-taxidlist`, use v5 with the required taxonomy data and inspect stderr and resulting taxids; the source reports an unfiltered successful exit when data are missing. See `references/database.md` for the caveat and `references/search.md` for filtering examples. Reciprocal best hits are only ortholog candidates; use `bio-ortholog-inference` for their limitations.

Done when requested sequences or filtered hits are checked against the supplied IDs and taxonomic scope.
