# Build a database

```bash
makeblastdb -in reference.fasta -dbtype nucl -blastdb_version 5 -parse_seqids -hash_index -out custom_db
```

Use `prot` for proteins and a fresh output prefix. `-parse_seqids` permits accession extraction; provide `-taxid` or `-taxid_map` when assignments are known. For prebuilt resources use `update_blastdb.pl --showall pretty` and check storage and release before download. v5 taxonomy filtering needs additional taxonomy files; see `references/database.md`.

Done when `blastdbcmd -db custom_db -info` confirms the database and its input and version are recorded.
