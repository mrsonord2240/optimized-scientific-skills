# Matched molecular pairs

```bash
mmpdb fragment data.smi -o data.fragments && mmpdb index data.fragments -o data.mmpdb && mmpdb loadprops -p props.tsv data.mmpdb
```

props.tsv needs ID values matching data.smi and numerical endpoint columns. Query and explicitly sort property-prefixed results as in `references/mmp.md`; interpret pair count, dependence and chemical context without a universal confidence cutoff.

Done when pair database, transformations and supporting evidence are reported.
