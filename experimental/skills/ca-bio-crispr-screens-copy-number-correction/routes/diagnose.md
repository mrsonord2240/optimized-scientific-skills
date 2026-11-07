# Detect CN bias
~~~bash
python scripts/detect_cn_bias.py gene_lfc.tsv copy_number.tsv diagnostic.json
~~~
TSVs contain gene, lfc and gene, copy_number respectively. Run genome-wide and for each amplicon with diploid background, before and after correction. The source uses CN >4 and diploid 1.5-2.5; rho <-0.1 with p<0.01 or gap <-0.5 with p<0.01 flags bias. Fewer than eight amplified genes can be underpowered: inspect gap even when nonsignificant. See [diagnostic details](../references/diagnostic.md).
Done when global and focal results, sample sizes and residual limitations are reviewed.
