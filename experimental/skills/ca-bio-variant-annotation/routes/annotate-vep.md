# Annotate with VEP on MANE Select

```bash
bash scripts/annotate_vep.sh input.vcf.gz reference.fa out/sample --cache-dir ~/.vep --fork 4
```

Normalizes (left-align, split multiallelics), then runs `vep --offline --assembly GRCh38 --everything --mane_select --pick --pick_order mane_select,mane_plus_clinical,canonical,biotype,rank`. Writes `out/sample_vep.vcf.gz` and a summary of VEP version, record count and HIGH/MODERATE count. It refuses any assembly but GRCh38 and an output with no `MANE_SELECT` field in CSQ (VEP 110 or later is needed).

- Add predictor plugins by repeating `--plugin`, one calibrated predictor per evidence type: `--plugin "dbNSFP,dbNSFP4.3a.gz,REVEL_score"`, `--plugin "AlphaMissense,file=AlphaMissense_hg38.tsv.gz"`, `--plugin "SpliceAI,snv=spliceai_snv.vcf.gz,indel=spliceai_indel.vcf.gz"`. Which to use: `routes/predictors.md`.
- The pick order is constrained so MANE leads and transcript length or accession never decides. Unconstrained `--pick` is clinically dangerous; see `routes/transcript-choice.md`.
- `--everything` turns on `--canonical`; report on the MANE block, not the canonical one.
- Flatten the result with `routes/parse-csq.md`.
- Record VEP version, cache version and plugin versions with the result.

Done when `out/sample_vep.vcf.gz` exists and the pinned build, transcript set, VEP version and cache version are stated.
