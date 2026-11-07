# RNA-seq length bias (GOseq)

```bash
Rscript scripts/goseq_ora.R de_results.csv results.csv gene_col=gene_id padj=0.05 lfc=1
```

Human hg38 only, Ensembl gene IDs (version suffixes stripped by the script). Gene lengths come from the local `TxDb.Hsapiens.UCSC.hg38.knownGene` and the GO map from the local `org.Hs.eg.db`; nothing is downloaded. The script stops if fewer than 10 genes map to both, and warns above 15% loss.

- Standard ORA on an RNA-seq list over-reports long-gene terms; RPKM or TMM does not fix that. GOseq (Wallenius) does.
- GOseq p-values are unadjusted: the script adds `padj` (BH). Use `padj`.
- Do not call `nullp(de_genes, 'hg38', 'ensGene')`; it can reach the network.
- Another organism or release: supply a positive named gene-length vector and a local GO-term-to-gene list in the same ID namespace, and record the TxDb and OrgDb versions.
- GSEA on a moderated t or Wald z statistic is largely immune to this bias (`routes/other-skills.md`).

Done when `results.csv` has `padj` and the mapping rates in `results.csv.log` are reported.
