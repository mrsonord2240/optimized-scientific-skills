# FOCUS gene fine-mapping
~~~sh
focus finemap gwas.sumstats 1000G_EUR_all focus.db --p-threshold 5e-8 --tissue Whole_Blood --locations 38:EUR --out gwas_focus
~~~
Provide underlying GWAS summary statistics, an ancestry-matched LD fileset covering analysed chromosomes, a matching model DB, and the correct genome-build/population locations. On Windows use relative positional paths because source-tested pyfocus splits colon-delimited arguments. The model DB can be built from a validated weight table:
~~~sh
python scripts/build_focus_db.py panel.tsv custom_focus.db --ref-name custom_panel --tissue Whole_Blood
~~~
Supply real per-gene cv_r2 and cv_r2_pval; source defaults are placeholders, not observed cross-validation. Assess prior sensitivity and report the full credible set; PIPs are model-dependent. MA-FOCUS requires matched inputs per ancestry and shared-effect assumptions. Source-reported version patches and LD path constraints are in [FOCUS notes](../references/focus-finemapping.md); verify installed software.
**Done when:** full locus results and credible sets are examined, including prior, panel coverage and ancestry assumptions.
