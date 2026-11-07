# Functional consequences

```bash
python scripts/hmmscan_to_pfamscan.py pfam_domtbl.txt pfam_scanfmt.txt Pfam-A.clans.tsv.gz
```

Generate the domain input with hmmscan --cut_ga and obtain matching Pfam clans first. Follow the ORF, sequence extraction and annotation order in `references/consequences.md`. Preserve annotated ORFs; request only imported consequence types, and use removeNoncodinORFs=FALSE for NMD. NMD predictions are probabilistic; see `references/nmd.md`.

Done when requested available consequences are annotated and their limitations stated.
