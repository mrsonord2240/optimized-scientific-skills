# Python CoMut plot

```bash
python scripts/comut_plot.py --mutations mutations.tsv --cohort cohort.tsv --clinical clinical.tsv --tmb tmb.tsv --output comut.pdf
```

Mutation, clinical and TMB TSVs use `sample`, `category`, `value`; the cohort TSV has `sample`. Clinical and TMB flags are optional. Supported alteration classes are Missense, Truncating, Splice, Amp, HomDel and Fusion. IDs are stripped and checked; TMB must be finite and nonnegative. The full cohort is retained; samples sort by TMB or displayed-mutation burden. Default top genes is 20 and sample labels hide above 50; use `--top` or `--sample-label-threshold` to change. The source reports comut 0.0.3 needs pandas 2.x; verify the installed contract. Output may be overwritten.

Done when class mapping, sample order, tracks, legend, cohort coverage and figure are inspected.
