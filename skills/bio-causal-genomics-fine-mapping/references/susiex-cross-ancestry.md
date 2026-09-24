# SuSiEx Cross-Ancestry Fine-Mapping

## Cross-Ancestry Fine-Mapping with SuSiEx

**Goal:** Jointly fine-map a locus across multiple ancestries assuming shared causal variants but population-specific LD.

**Approach:** Per-ancestry summary statistics + per-ancestry PLINK reference panel; SuSiEx computes each population's in-sample LD itself (from `--ref_file`, writing it to the `--ld_file` prefix — this is an output path, not a pre-built matrix you supply) and runs a joint SuSiE model with the resulting population-specific R matrices. SuSiEx assigns populations by the ORDER of the comma-separated `--sst_file`/`--n_gwas`/`--ref_file`/`--ld_file` lists (there is no `--pop` flag); keep all four lists in the same population order.

Full command with every required column flag (`--chr_col` ... `--pval_col`, one value per population), `--ld_file`, `--plink`, `--level` and output layout: `examples/susiex_multiancestry.sh` (also the one-line form in `SKILL.md`).

**Checked on SuSiEx 1.1.2** (`SuSiEx --help`): `--plink=<path/to/plink binary>` is required — omitting it is not caught by argument parsing but fails later when SuSiEx shells out to PLINK to build each population's LD matrix. Verified end-to-end on a synthetic two-population locus with a planted shared causal SNP: SuSiEx recovered it at CS_PIP = 1.0 in a single credible set.

The output includes per-population PIPs and a joint credible set (`<out_name>.cs`, `<out_name>.snp`, `<out_name>.summary`). Credible sets from SuSiEx are typically 2-5x smaller than EUR-only susie_rss when AFR is included, because AFR shorter LD blocks resolve EUR-tagged regions.
