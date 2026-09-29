# RASQUAL joint modeling

Use this reference only for the moderate-cohort branch that combines total
accessibility with allele-specific counts.

RASQUAL uses a non-standard command line. It reads VCF records from standard
input and uses single-letter flags. It does not accept the `--features`,
`--counts`, or `--vcf` flags used by some newer wrappers.

## Prepare inputs

Use `rasqualTools::saveRasqualMatrices` from R to create the binary count and
offset files consumed by `-y` and `-k`. Prepare one row per feature with:

1. feature name;
2. chromosome;
3. start;
4. end;
5. number of testing SNPs in the cis window; and
6. number of feature SNPs overlapping the peak.

Derive the two SNP counts from the exact VCF records that will be streamed for
that feature. Skip a feature with zero feature SNPs.

## Invoke one feature at a time

```bash
python3 scripts/run_rasqual_features.py \
  --features features.tsv \
  --counts counts.bin \
  --offsets offsets.bin \
  --vcf cohort.vcf.gz \
  --samples "$N_SAMPLES" \
  --matrix-rows "$N_FEATURE_ROWS" \
  --output-dir rasqual_results
```

The driver validates the six-column feature table, requires each binary matrix
to contain exactly `samples * matrix_rows` doubles, checks the VCF header sample
count, assigns each feature line's 1-based matrix row to `-j`, skips and records
zero-feature-SNP rows, and refuses to reuse an output directory. Use `--dry-run`
to inspect the exact argument arrays without executing the tools.

Flag meanings in the provider workflow:

- `-y`: binary count file.
- `-k`: binary size-factor or offset file.
- `-n`: number of samples.
- `-j`: row index of the feature.
- `-l`: number of testing SNPs in the cis-window.
- `-m`: number of feature SNPs overlapping the feature.
- `-s`, `-e`, `-f`: feature start, end, and name.

RASQUAL estimates genotype and allelic correlation internally from the
tabix-streamed VCF. Do not add an external LD-precomputation step.

Inspect `rasqual --help` for the installed checkout and prefer the current
`rasqualTools` wrapper where it is available. Record the checkout commit and
matrix-generation procedure because the interface is unusual.

The pinned RASQUAL output exposes a likelihood-ratio chi-square in column 11.
The driver converts each reported one-degree-of-freedom feature-variant test to
a p-value and applies BH across the complete reported test family in
`association_summary.tsv`. If another checkout changes that schema, pass a
verified `--chi-square-column`; do not silently reuse the default. Power depends
on the study design, depth, allele frequency, and model fit, so this Skill makes
no universal fold-gain claim.

