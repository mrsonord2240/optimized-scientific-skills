# hashsolo
~~~sh
python scripts/hashsolo_run.py gex.h5ad hto_counts.csv HTO_A,HTO_B,HTO_C,HTO_D demuxed.h5ad
~~~
The HTO CSV has cell barcodes as its first column and raw tag-count columns; all query barcodes must match. Source priors are negative/singlet/doublet = 0.01/0.8/0.19 and noise barcodes are set to 1 for two or three tags. These are assumptions, not universal loading rates. Check high Negative fractions, tag balance and actual API behavior; retain all classifications.
**Done when:** calls, priors, tag counts and unassigned/doublet fractions are reviewed.
