# Add to an alignment

```bash
mafft --add new_sequences.fasta existing_alignment.fasta > updated.fasta
```

Use `--addfragments` for partial reads, `--addprofile` for an aligned profile, and `--keeplength` only when insertions must be dropped to preserve existing columns. ClustalOmega and HMMER alternatives are in `references/profiles.md`.

Done when the updated alignment is written and expected IDs and coordinate changes are checked.
