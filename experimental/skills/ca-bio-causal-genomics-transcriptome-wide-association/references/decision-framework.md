## TWAS - MR - Coloc Triangulation

A TWAS-significant gene is associational, not causal. The strongest defensible claim that a gene mediates a GWAS effect comes from triangulating three orthogonal lines of evidence:

| Line | What it tests | Threshold |
|------|---------------|-----------|
| TWAS | Predicted-expression association with trait | S-MultiXcan joint p < 2.3e-6, OR S-PrediXcan per-tissue p < 4.6e-8 (49-tissue Bonferroni), OR per-tissue FDR < 0.05 |
| cis-eQTL MR | Causal effect of expression on trait under IV assumptions (cross-reference bio-causal-genomics-mendelian-randomization) | Wald-ratio or IVW p < 0.05/n_genes; instrument F > 10 |
| Colocalization | Shared causal variant between GWAS and eQTL (cross-reference bio-causal-genomics-colocalization-analysis) | coloc.abf or coloc.susie PP.H4 >= 0.7 |
| FOCUS | Probabilistic per-gene PIP under TWAS fine-mapping | PIP >= 0.8 |

**Operational rule:** Report a gene as a "strong candidate causal gene" only when 3 of the 4 are concordant (TWAS hit + coloc PP.H4 >= 0.7 + FOCUS PIP >= 0.8, with cis-MR as a supporting fourth). 2-of-4 concordance is "suggestive"; 1-of-4 is "associational only". The combination is more conservative than any single method but matches the standards used in modern GWAS-to-target pipelines (Open Targets Genetics Mountjoy 2021 Nat Genet 53:1527; FinnGen R10 release notes).

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|--------------------|
| S-MultiXcan joint p (gene-wide) | < 2.3e-6 (0.05 / 22k genes) | One p per gene across tissues; standard genome-wide TWAS significance |
| S-PrediXcan per-tissue (49 tissues) | < 4.6e-8 (0.05 / (49 x 22k)) | Cross-tissue x cross-gene Bonferroni |
| Per-tissue FDR (alternative) | BH q < 0.05 within each tissue | Less conservative; report cross-tissue replication of FDR-significant hits |
| FUSION cross-validation R^2 | >= 0.01 (per gene) | FUSION default heuristic; below this, gene is not heritable in tissue and weights drop |
| FUSION heritability p | < 0.01 (hsq_p) | FUSION default; filters out non-heritable expression |
| Coloc PP.H4 for triangulation | >= 0.7 | Open Targets / common practice; >= 0.8 for high-confidence |
| FOCUS gene PIP (causal) | >= 0.8 | Mancuso 2019 convention; PIP >= 0.5 is suggestive |
| cis-eQTL MR instrument F | >= 10 | Standard MR convention to avoid weak-instrument bias |
| Tissue eQTL N (weight stability) | >= 100 donors | Below this, per-gene elastic-net weights unstable |
| cis-window radius | +/- 500 kb of gene TSS/TES (FUSION) or +/- 1 Mb (S-PrediXcan) | Captures most cis-eQTL signal; window choice rarely changes top hits |
| LD reference panel size | >= 500 individuals matched ancestry | 1000 Genomes superpopulation reference is standard |
| MA-FOCUS minimum ancestry N | >= 2 ancestries with significant TWAS | Joint inference requires non-trivial heterogeneity |
| S-MultiXcan condition number | < 30 (correlation matrix) | PCA regularisation kicks in above this; near-collinear tissues collapsed |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Error: weight file <gene>.RDat not found` (FUSION) | weights.pos points to relative paths; `--weights_dir` not set correctly | Verify `--weights_dir` matches the .pos `WGT` column directory |
| `KeyError: 'SNP'` or `'A1'` (S-PrediXcan) | GWAS sumstat columns not mapped via `--snp_column` etc. | Pass `--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA --pvalue_column P` explicitly |
| All TWAS Z's are 0 or NA | Allele coding mismatch between GWAS and weight panel | Harmonise (effect-allele flip); check `--keep_non_rsid` if using non-rsID SNPs |
| Genomic inflation lambda >> 1.05 | Wrong LD reference OR population structure not controlled in upstream GWAS | Fix at the GWAS stage; do not adjust TWAS lambda post-hoc |
| FOCUS reports PIP = 1 for one gene at every locus | Only one gene in panel at locus; degenerate posterior | Expand panel coverage or report locus as panel-limited |
| S-MultiXcan condition number warning | Tissues near-collinear (multiple brain regions) | Increase regularisation (`--regularization 0.5`) or restrict to tissue subset |
| FOCUS DB not matching FUSION weights | Custom weight panel without FOCUS DB | Build one (`focus import`, or the direct-build script `scripts/build_focus_db.py`, see `references/focus-finemapping.md`) |
| `DataFrame.pivot() takes 1 positional argument but 4 were given` or `'DataFrame' object has no attribute 'append'` (`focus finemap`, at "Calculating PIPs") | pyfocus 0.802 `finemap.py:1012` / `:188` use APIs removed in pandas 2.0 | Apply the two `finemap.py` patches in version-notes.md |
| `focus import` finishes with an empty DB / 0 genes imported | `mygene` or `rpy2` missing (only an ERROR log line); rpy2 also needs R built as a shared library | Install both, or build the DB directly (`scripts/build_focus_db.py`) |
| MA-FOCUS H0 probability dominates | Cross-ancestry heterogeneity at the locus | Run per-ancestry FOCUS separately; do not force joint |
| S-PrediXcan output has effect sizes much larger than expected | `sdY` proxy mis-specified; standardised vs unstandardised mismatch | Confirm GWAS Z and beta scale; rerun with `--additional_output` |
| `InvalidArguments: Specify either cutoff_ratio or cutoff_threshold` (S-MultiXcan) | `--cutoff_condition_number` (or `--cutoff_ratio`/`--cutoff_threshold`) omitted; not optional in installed versions | Always pass `--cutoff_condition_number 30` |
| A comma-delimited filter over S-MultiXcan's `--output *.csv` silently returns 0 rows (no error) | `SMulTiXcan.py` hard-codes tab-separated output regardless of the output filename's extension | Filter with `awk -F'\t'` (see `routes/metaxcan.md`), not `-F','` |
| `TypeError: read_csv() got an unexpected keyword argument 'delim_whitespace'` (`focus`/`focus finemap`) | Unpinned pandas>=2.2 | Install pins and patch in version-notes.md |
| `AttributeError: module 'numpy' has no attribute 'warnings'` (`focus`/`focus finemap`) | Removed `np.warnings` alias; pinning pandas does not help | Patch after install (version-notes.md) |
| `Please specify independent regions location or default regions with '37:EUR', etc.` (`focus finemap`) | `--locations` omitted | Always pass it (see `references/focus-finemapping.md`) |
| `focus finemap` reports an implausible population count (e.g. "Detecting 2 populations" from one file) on Windows | Absolute Windows paths get colon-split by pyfocus's multi-ancestry parser | Use relative paths (see `references/focus-finemapping.md`) |
| `FileNotFoundError: ...bim` on the `ref` argument to `focus finemap` | LD reference path templated like FUSION's per-chromosome `--ref_ld_chr` (e.g. `.../chr` with no matching file); pyfocus passes it straight to `pandas_plink.read_plink()`, which needs an exact prefix or a real glob, not a chr-substitution convention | Point at one PLINK fileset covering every chromosome analyzed, or a `pandas_plink` glob (`prefix/chr*.bed`); `--chr`/`--locations` filter after loading, they do not select the file (see `references/focus-finemapping.md`) |
| A PIP filter against a `pip` column silently returns nothing (or the whole row) | Installed pyfocus 0.802 never writes a plain `pip` column -- output is `pips_pop1` (population-indexed even for a single population) or `pips_me` for the cross-ancestry marginal PIP | Filter on `pips_pop1` (single-ancestry) or `pips_me` (MA-FOCUS); see `routes/focus.md` |
| `Error in wgt.matrix[qc$flip, ] : incorrect number of dimensions` or `non-conformable arguments` (`FUSION.post_process.R`) | Single-SNP ("top1") weight models: `wgt.matrix[m.keep,]` and `genos$bed[,m[m.keep]]` (lines ~168/170 and ~251/254) drop a 1-row/1-column matrix to a bare vector without `drop=FALSE` | Patch your local `fusion_twas` clone: add `,drop=FALSE` to both subsetting lines at both locations, or filter single-SNP genes out of `--input` before running `post_process.R` |
