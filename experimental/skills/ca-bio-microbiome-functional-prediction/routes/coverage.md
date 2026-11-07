# Report NSTI coverage

```bash
python scripts/nsti_report.py picrust2_out asv_table.tsv --max-nsti 2.0
```

The TSV has ASVs in rows and samples in columns, with no leading BIOM comment; the script expects the source-version combined NSTI file and prints mean, median and the ASV/read fraction above the threshold among IDs present in that file. It does not count input ASVs absent from the NSTI file, failed placements or missing NSTI values as dropped; separately reconcile the original input, placement results and missing values before reporting total coverage. Reference distance is an uncertainty indicator, not validation of predicted function.

Done when the NSTI distribution and all exclusion categories and affected read fractions are reconciled and reported.
