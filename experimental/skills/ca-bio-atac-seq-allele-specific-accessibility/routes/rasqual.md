# Joint cohort model

```bash
python scripts/run_rasqual_features.py --features features.tsv --counts counts.bin --offsets offsets.bin --vcf cohort.vcf.gz --samples 50 --matrix-rows 1000 --output-dir rasqual_results
```

Use WASP-corrected inputs. Prepare binary matrices and six-column features with testing and feature SNP counts from the actual cis window, as in `references/rasqual-joint-modeling.md`. Verify installed RASQUAL output schema; column 11 is the source default. The driver reports BH across the complete reported feature-variant family. Set samples and rows to actual inputs.

Done when joint results and defined multiplicity are reviewed.
