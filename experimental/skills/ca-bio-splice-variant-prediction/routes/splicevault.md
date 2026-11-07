# Empirical splice outcomes

```bash
python scripts/splicevault_lookup.py chr17 7674292 T C --tsv local.tsv.gz --transcript ENST00000269305 > events.txt
```

Use a local indexed GRCh38 table to keep coordinates local. Remote lookup sends coordinates to a third party; use it only when appropriate. Empty output is not benign evidence. See `references/splicevault.md` for event interpretation.

Done when matching transcript events or the absence of a matching record are reported.
