# Diagnose copy-number bias
~~~bash
python scripts/cn_bias.py gene_lfc.tsv copy_number.tsv > cn_bias.txt
~~~
Provide gene, lfc and gene, copy_number TSVs. Source tests genome-wide rho and amplified-versus-diploid gap (CN >4 versus 1.5-2.5), with p<0.01 and gap <-0.5 for focal bias. Repeat for a candidate region with diploid background, and interpret limited sample power: a nonsignificant result does not establish absence of bias. See [CN notes](../references/copy-number-bias.md) and canonical bio-crispr-screens-copy-number-correction.
Done when global and focal results, counts and uncertainty are reviewed.
