# Compare predictor outputs

```python
from scripts.splice_parsers import build_concordance
T = build_concordance('input.vcf', spliceai_vcf='output.vcf', pangolin_vcf='pangolin_output.vcf', mmsplice_csv='mmsplice_predictions.csv')
T.to_csv('concordance.tsv', sep='\t')
```

The parser is `scripts/splice_parsers.py`. Split multiallelic inputs before Pangolin; retain and report n_scored and missing outputs. Agreement is prioritization, not proof. For additional predictors see `references/other-predictors.md`; interpretation and limitations are in `references/interpretation.md`.

Done when per-variant scores, missingness and disagreement are reported.
