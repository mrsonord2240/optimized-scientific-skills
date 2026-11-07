# Sequence records

```bash
python scripts/fetch.py --email researcher@example.org --db nucleotide --id NM_007294.4 --rettype gb --retmode text --out record.gb
```

Choose `fasta` for sequence, `gb` for nucleotide features, `gp` for protein features, or `gbwithparts` for assembled WGS content when appropriate. `fasta_cds_aa` retrieves translated CDS records. Verify response and requested accession.version; empty or HTML error bodies are not records. For large lists use `bio-batch-downloads`. Other format choices and parsing notes are in `references/formats.md`.

Done when the requested record is saved and its accession, format and content checked.
