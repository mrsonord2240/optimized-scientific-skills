# Splice-site sequence strength

```bash
python scripts/splicing_qc.py sites --donor CAGGTAAGT --acceptor TTTTTTTTTTTTTTCCTTAGGAG > sites.txt
```

Supply 9-nt donor or 23-nt acceptor sequences; invalid sequences produce NaN. MaxEnt and contextual SpliceAI answer different questions; see `references/sites.md`.

Done when requested sequence scores and limitations are reported.
