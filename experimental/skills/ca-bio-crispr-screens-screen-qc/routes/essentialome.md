# Evaluate reference recovery
~~~bash
python scripts/essentialome_recovery.py gene_lfc.tsv CEGv2.txt NEGv1.txt > essentialome.txt
~~~
Gene LFC TSV has gene and lfc columns; reference sets must match the assay context and contain both positive and negative genes. The helper reports PR-AUC and ROC-AUC based on -LFC; its >0.7 and <0.5 labels are source conventions, not universal pass/fail proofs. Review overlap, class prevalence and biological context. See [reference recovery](../references/essentialome-recovery.md).
Done when reference coverage and PR/ROC results are interpreted in context.
