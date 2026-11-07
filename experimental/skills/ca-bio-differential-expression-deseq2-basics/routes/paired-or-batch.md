# Batch, or the same subjects in both groups

```bash
Rscript scripts/deseq2_de.R counts.csv samples.csv results.csv ref=control block=batch
```

`block=` names the nuisance column in `samples.csv`: the batch, or the subject ID when each subject was measured in both groups. Several columns: `block=batch,sex`.

- Put the batch or subject in the model with `block=`. Never subtract it from the counts first.
- Paired samples: `block=<subject column>`. Each subject then acts as its own control.
- `design matrix not full rank` means a block level sits entirely inside one group, so the two cannot be separated. Drop that covariate or report that the effect is confounded.
- Significant genes are the rows with `padj < 0.05` in `results.csv`.

Done when `results.csv` exists and you have answered from it.
