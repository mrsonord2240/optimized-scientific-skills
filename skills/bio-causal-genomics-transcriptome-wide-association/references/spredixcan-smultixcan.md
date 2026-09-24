# S-PrediXcan + S-MultiXcan Pipeline

**Goal:** Run TWAS across all GTEx tissues using pre-trained PredictDB models, then combine via S-MultiXcan for a joint multi-tissue test.

**Approach:** For each tissue, run S-PrediXcan with the matched model DB and covariance file; collect per-tissue outputs into a folder; run S-MultiXcan with the same model folder and GWAS to produce a joint multi-tissue Z and per-tissue significance.

PredictDB models are at predictdb.org; GTEx v8 MASHR-EUR is the standard EUR panel, with a `.db` (model) and `.txt.gz` (covariance) file pair per tissue. The runnable pipeline is `examples/s_predixcan_pipeline.sh`: per-tissue `SPrediXcan.py` (`--model_db_path`, `--covariance`, `--gwas_file`, explicit `--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA --pvalue_column P`), then `SMulTiXcan.py` (`--models_folder`, `--models_name_pattern`, `--snp_covariance`, `--metaxcan_folder`, `--metaxcan_filter`, `--metaxcan_file_name_parse_pattern`, `--cutoff_condition_number 30`), then a filter on the joint p at 2.3e-6.

`--cutoff_condition_number 30` (the canonical MetaXcan setting) is required, not optional; omitting it raises the `InvalidArguments` error in SKILL.md Common Errors (confirmed 2026-09-19 against a real run).

S-MultiXcan applies PCA regularisation on the inter-tissue correlation matrix: `--cutoff_condition_number 30` drops near-collinear components, and `--regularization 0.1` (off unless passed explicitly) adds a ridge. Tissues that are nearly collinear with another (e.g. multiple brain sub-regions) are absorbed into shared components and do not contribute independent power.
