# One factor with three or more levels

One level against another:

```bash
Rscript scripts/deseq2_de.R counts.csv samples.csv results.csv condition=dose ref=none versus=high
```

Does the factor change the gene at all, across every level:

```bash
Rscript scripts/deseq2_de.R counts.csv samples.csv results.csv condition=dose ref=none test=lrt
```

- Run the first form once per pair you need. Each run writes its own file.
- The `test=lrt` output has no fold-change column. It answers "any difference", so take effect sizes from a `versus=` run.
- Add `block=<column>` for a batch or subject, as in `routes/paired-or-batch.md`.

Done when each requested comparison has its results file.

More detail on what the LRT does and does not report: `references/lrt.md`.
