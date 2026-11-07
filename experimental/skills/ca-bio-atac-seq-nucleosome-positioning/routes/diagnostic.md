# Fragment and positioning diagnostic

```bash
python scripts/vplot.py sample.bam tss.bed vplot.png
```

The BED start is the feature center; supply strand in column 6 to mirror minus-strand features. If the aggregate is flat, inspect annotation and library quality before calling positions. `python scripts/estimate_nrl.py sample.bam` reports the 150–250 bp mono-fragment mode, an approximate repeat length rather than a precise measurement. The source windows are NFR <100, mono 180–247 and di 315–473 bp; adjust to the organism. For hg38 R-based class summaries and TSS heatmap use `Rscript scripts/nucleosome_analysis.R sample.bam prefix tss.bed`. See `references/interpretation.md`.

Done when the diagnostic is inspected and fragment windows, features and limitations are reported.
