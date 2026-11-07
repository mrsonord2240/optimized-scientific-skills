# Annotation-free intron usage

```bash
Rscript scripts/diff_splicing_leafcutter.R
```

First edit sample names and settings at the top; set LEAFCUTTER_DIR to the source clone and run where BAMs are available. BAMs require XS tags. Nonstandard contigs require NONSTANDARD_CONTIGS=TRUE; lower cluster minimum for shallow data. For two or three replicates treat calls as exploratory and seek concordance. Batch covariates and cluster interpretation are in `references/leafcutter.md` and `references/interpretation.md`.

Done when cluster significance and intron effects are available and interpreted.
