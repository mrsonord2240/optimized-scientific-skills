# Signac processing

```bash
Rscript scripts/signac_workflow.R filtered_peak_bc_matrix.h5 fragments.tsv.gz singlecell.csv 4 1000
```

Arguments four and five are minimum TSS enrichment and fragment count. Use ARC `per_barcode_metrics.csv` for matched Multiome input; the script selects the Peaks matrix and uses GEX-style barcodes. It assumes hg38 and writes `scatac_signac.rds` and diagnostic PDFs in the current directory; the CLI does not expose an output prefix. It does not run doublet detection. Defaults and the ARC blacklist-ratio caveat are in `references/qc.md`; adjust thresholds deliberately for shallow data. LSI dimensions 2:30 omit the depth-associated first component.

Done when QC, retained cells, embedding and gene activity are checked and assumptions and filter changes are reported.
