# Score variants

```bash
python scripts/score_variants_torch.py chrombpnet_nobias.h5 genome.fa variants.tsv scores.tsv
```

Require passing model QC first. Input is five headerless tab-separated fields: chromosome, 1-based position, REF, ALT, ID. Script uses the forward strand and log-count head, reporting log2 count change. It does not produce null p-values. For null-based p-values and reverse-complement averaging use variant-scorer as described in `references/method-reference.md`; correct across variants and report effect and model quality. A single-tool call is exploratory; compare appropriate orthogonal predictions.

Done when scored variants, skipped counts, model QC and interpretation are reported.
