# Build or reuse a gzip FASTA index
```bash
python scripts/pyfastx_index.py reference.fasta.gz --record-id seq_00042
```
**Done when:** JSON reports the requested ID, count and length and the persistent .fxi sidecar exists. The optional ID defaults to the first indexed record. Confirm the index is appropriate for the current input before reuse.
For multi-file Biopython `SeqIO.index_db`, BGZF and plain input are supported but plain gzip is not; indexed IDs must be unique. For linear FASTQ and quality encodings, and pyfastx alternatives, see [reader selection](../references/readers-and-indexing.md).
