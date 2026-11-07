# Known identifiers

```bash
epost -db nucleotide -input ids.txt | efetch -format fasta > output.fasta
```

Supply one valid ID per line and verify resulting records against the requested list. For large lists, source EPost chunking and per-key cursor details are in `references/ids.md`. Use `-format fasta`; `-mode webenv` is not a history-server selector. Avoid simultaneous streams without a shared request-rate limiter.

Done when the requested identifiers are recovered or missing IDs are explicitly reported.
