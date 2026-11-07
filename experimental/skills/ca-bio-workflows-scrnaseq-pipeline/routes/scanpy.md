# Single-sample Scanpy
~~~bash
python scripts/scanpy_workflow.py filtered_feature_bc_matrix.h5 results_scanpy
~~~
Done when QC, doublet calls, clustering and marker rankings are reviewed. Output is a new directory; the object is clustered, not biologically annotated. This illustrative filtered-input workflow does not remove ambient RNA, demultiplex, integrate samples or run sample-level DE. Apply those handoffs separately. Flat 200/5000 gene and 20% mitochondrial limits are examples; adapt to tissue and per-sample distributions. Scanpy gene-expression-only import excludes other modalities. Raw count layer is subset to HVGs in the saved object; retain the original full counts for pseudobulk.
