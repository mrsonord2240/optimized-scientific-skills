# GO-slim breakdown

```bash
Rscript scripts/groupgo.R genes.txt results.csv keytype=ENTREZID ont=BP level=3
```

`groupGO` classifies genes at one fixed DAG level. It is not a test: the output is counts per term and has no p-value. Never read a count as significance, and use `routes/ora.md` when the question is what is over-represented.

Done when `results.csv` exists and the answer says it is descriptive.
