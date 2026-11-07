# Independent site concordance

```bash
bash scripts/site_concordance.sh bound.bed unbound.bed other_calls.bed regions.bed
```

Omit the fourth argument if no common-region restriction is needed. For independent calls, HINT-ATAC and Wellington commands are in `references/inputs.md`. Compare overlap fractions at bound and unbound sites and the reverse overlap; no universal concordance threshold is established. Single-tool results remain exploratory.

Done when bound and unbound overlap and reverse overlap are reported with the region restriction and method.
