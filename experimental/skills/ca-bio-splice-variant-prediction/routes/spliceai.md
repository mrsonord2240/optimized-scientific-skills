# SpliceAI prediction

```bash
python scripts/spliceai_clingen_classify.py input.vcf genome.fa --build grch38 --out-prefix results
```

Use grch37 when appropriate. The reference is human and must match the VCF; check REF before running. Default D=50 and M=0; scored variants below 0.20 are rerun at D=500. Missing or unsupported scores remain not_scored. The labels are supporting computational evidence only, and gene-level scores do not select a disease-relevant transcript. See `references/spliceai.md` and `references/selection.md`.

Done when scored and unscored variants are reported with version, reference, window and transcript context.
