# S-PrediXcan
~~~sh
python SPrediXcan.py --model_db_path tissue.db --covariance tissue.txt.gz --gwas_file gwas.tsv --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA --pvalue_column P --output_file tissue.csv
~~~
Use matched tissue model and covariance, harmonized alleles and ancestry; execute in the MetaXcan software environment. To combine a prepared set of tissue outputs, use the source's S-MultiXcan settings:
~~~sh
python SMulTiXcan.py --models_folder models --models_name_pattern 'mashr_(.*)\.db' --snp_covariance covariance.txt.gz --metaxcan_folder results --metaxcan_filter '(.*)\.csv' --metaxcan_file_name_parse_pattern '(.*)\.csv' --gwas_file gwas.tsv --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA --pvalue_column P --regularization 0.1 --cutoff_condition_number 30 --verbosity 9 --throw --output joint.tsv
~~~
Select tissues in advance and use one consistent model family; the source's 14-tissue example is representative, not all 49. S-MultiXcan output is tab-separated in the source-tested version even when named .csv. Verify the installed interface. [Methods](../references/spredixcan-smultixcan.md), [tissue choice](../references/tissue-and-model-selection.md).
**Done when:** per-tissue and, if requested, joint output passes harmonization and multiple-testing review.
