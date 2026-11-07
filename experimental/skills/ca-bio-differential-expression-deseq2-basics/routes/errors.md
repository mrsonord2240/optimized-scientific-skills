# Errors and results that look wrong

| Symptom | Cause | Fix |
|---------|-------|-----|
| `design matrix not full rank` | A covariate is confounded with the condition or another covariate | Drop it, or report the effect as inseparable |
| `counts are not integers` | Salmon/kallisto/RSEM estimates, or normalized values | `routes/tximport.md`; normalized values cannot be used |
| `sample IDs differ` | Count columns and the sample table disagree | Make the first column of `samples.csv` match the count column names |
| Fold changes have the opposite sign to expected | Wrong `ref=` | Rerun with the right reference level |
| Nearly every gene is significant in one direction | Most of the transcriptome really shifted, which breaks the default normalization | `references/size-factors.md` |
| A dispersion plot shows the fitted line missing the points | The default trend does not fit this data | `references/dispersion-trend.md` |
| An old tutorial's `betaPrior` or `contrast=` code gives different numbers | Behaviour changed across DESeq2 versions | `references/betaprior-history.md` |
