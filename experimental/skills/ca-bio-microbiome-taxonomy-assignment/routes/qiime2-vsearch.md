# QIIME2 alignment-consensus (classify-consensus-vsearch)

```bash
bash scripts/qiime2_vsearch.sh rep-seqs.qza ref-seqs-region.qza silva-138-99-tax.qza taxonomy_out
```

Aligns each query globally with VSEARCH, keeps the top `MAXACCEPTS` (10) hits at or above `PERC_IDENTITY` (0.8), and reports a rank only if at least `MIN_CONSENSUS` (0.51) of them agree; otherwise that rank is dropped. `THREADS` defaults to 8. Writes `taxonomy-vsearch.qza` and `vsearch-hits.qza`.

- Stores no pickled model, so it is immune to the scikit-learn version error. Hits are transparent; it is slower than naive Bayes.
- Use a region-extracted reference for amplicon queries (`extract-reads` as in `routes/qiime2-region-nb.md` step 1).

Done when `taxonomy-vsearch.qza` exists and the identity and consensus values are recorded.
