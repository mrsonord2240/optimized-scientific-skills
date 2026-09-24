# FUSION Pipeline

**Goal:** Run sumstat TWAS using FUSION and conditional joint analysis to identify independently associated genes.

**Approach:** Format GWAS sumstats to FUSION expected columns (SNP A1 A2 Z); run `FUSION.assoc_test.R` per chromosome with pre-computed weights and ancestry-matched LD; post-process with `FUSION.post_process.R` (conditional/joint analysis run by default) to identify independent genes; flag conditional-significant genes for follow-up.

Pre-computed FUSION weights (e.g. GTEx v8 Whole_Blood) live at http://gusevlab.org/projects/fusion/: download the `.pos` summary and the per-gene RData files into the weights directory. GWAS sumstats need the columns `SNP A1 A2 Z` (Z on the standardised scale); use TwoSampleMR or a custom munger to harmonise alleles upstream.

Run `scripts/fusion_twas.sh` (in the Skill folder) from inside the `fusion_twas` clone (it loops `FUSION.assoc_test.R` over chromosomes, concatenates the tables, then runs `FUSION.post_process.R --plot --locus_win 100000` per chromosome; `twas_joint*.dat` reports per-gene conditional Z, and genes with joint Z > 4 are independent):

```bash
<skill-dir>/scripts/fusion_twas.sh gwas.sumstats gtex_whole_blood.pos gtex_whole_blood_wgt/ 1000G_EUR_LD/EUR. twas_out/twas "$(seq -s' ' 1 22)"
```

`FUSION.post_process.R` returns conditionally independent gene signals, not PIPs; use FOCUS for probabilistic fine-mapping. The `--locus_win 100000` parameter defines the conditioning window; 100 kb is conservative for non-HLA loci. FUSION's `--coloc_P` flag runs single-SNP coloc internally but is less robust than running coloc separately on the per-gene top eQTL.

**Known upstream crash on single-SNP ("top1") genes** -- see Common Errors in SKILL.md for the exact error, cause, and a verified two-line patch.
