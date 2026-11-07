# Import alignment output

```bash
python scripts/import_alignment.py AlignResult.mdalign output-prefix
```

The parser expects the source-described four metadata rows followed by a tab-separated header, with exactly one `Class` anchor in row one and sample columns after it. It writes `.features.tsv` and `.intensity.tsv`, retaining all features and sample values; existing outputs are overwritten. Confirm sample names and numerical conversion against the actual MS-DIAL version. It does not apply a Fill% filter or assign confidence. The source's common Fill% screen is 0.70 on its observed 0–1 scale, but verify that scale and group-wise missingness in the actual export. MS/MS assignment is not confirmation; apply further QC and evidence review with `bio-metabolomics-normalization-qc` and `bio-metabolomics-metabolite-annotation`.

Done when metadata and sample columns, identifiers, numeric values, missingness and output provenance are checked.
