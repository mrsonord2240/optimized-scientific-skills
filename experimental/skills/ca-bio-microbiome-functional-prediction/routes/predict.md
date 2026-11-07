# Predict functional potential

```bash
picrust2_pipeline.py -s asv_seqs.fna -i asv_table.biom -o picrust2_out -p 8 --hsp_method mp --max_nsti 2 --verbose
```

Supply matched representative-sequence and abundance inputs, a new output directory and a verified installed PICRUSt2/reference version. The source notes BIOM-converted TSV exports may contain a leading comment line that the pipeline cannot parse; use compatible input. Its example `examples/run_picrust2.sh` uses fixed file paths and runs prediction, descriptions and an NSTI summary, and requires adaptation. Predicted pathway abundance, optional coverage and per-sequence contributions are different outputs; `--per_sequence_contrib` requires `--stratified` in the source version. Compare placement and reference exclusions with the original table, then use [coverage](coverage.md).

Done when input provenance, version, placement, exclusions, NSTI and the prediction-only interpretation are recorded.
