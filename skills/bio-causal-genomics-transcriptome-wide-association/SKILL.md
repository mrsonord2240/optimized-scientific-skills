---
name: bio-causal-genomics-transcriptome-wide-association
description: Performs gene-level association from GWAS summary statistics via genetically predicted tissue expression using FUSION, PrediXcan, S-PrediXcan, S-MultiXcan, UTMOST, MOSTWAS, kTWAS, EpiXcan, TIGAR-V2, and probabilistic fine-mapping with FOCUS and MA-FOCUS. Use when running TWAS from GWAS sumstats, prioritising candidate causal genes from a GWAS lead locus, picking single-tissue vs cross-tissue models, identifying LD-induced TWAS false positives, choosing ancestry-matched prediction weights, fine-mapping co-regulated TWAS hits, or triangulating TWAS with cis-eQTL Mendelian randomization and colocalization to nominate a causal gene.
tool_type: mixed
primary_tool: FUSION
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: FUSION (head of `gusevlab/fusion_twas`, scripts dated 2023+), MetaXcan / S-PrediXcan / S-MultiXcan 0.7.5+ (`hakyimlab/MetaXcan`), PrediXcan model files from PredictDB (GTEx v8 elastic-net + MASHR), UTMOST (head of `Joker-Jerome/UTMOST`), pyfocus 0.8+ (`bogdanlab/focus`), MA-FOCUS (head of `mancusolab/ma-focus`), TIGAR-V2 (head of `yanglab-emory/TIGAR`), PLINK 1.9 + PLINK 2.0, R 4.3+, Python 3.9-3.11.

Before using code patterns, verify installed versions match. If versions differ:
- R: `Rscript --version`; for FUSION scripts inspect `--help` flags directly in the source
- Python: `pip show pyfocus` (MetaXcan is git-cloned, not on PyPI) then `SPrediXcan.py --help`, `SMulTiXcan.py --help`, `focus finemap --help`
- CLI: `plink2 --version`; FUSION ships as R scripts not a binary

If a script throws an error about an argument that has moved (e.g. `--gwas_file` vs `--gwas-file`) or a model database schema change, introspect the installed script with `--help` and adapt rather than retrying. PredictDB model file paths change with GTEx version; pin the version explicitly in scripts.

# Transcriptome-Wide Association

**"Find genes whose predicted tissue expression is associated with my GWAS trait"** -> Train SNP -> expression prediction models on a reference eQTL panel, apply the per-gene SNP weights to GWAS summary statistics or genotypes, and produce a gene-level Z-score equivalent to a weighted sum of SNP Z-scores. The output is a gene-by-tissue association, but TWAS is NOT direct evidence of causal mediation: an LD-tagged eQTL signal produces the same statistical association as a truly causal one, and the dominant failure modes are LD-induced false positives at gene-dense loci, tissue mis-specification, and ancestry mismatch between GWAS and prediction weights.

- CLI (sumstat TWAS, R): `FUSION.assoc_test.R --sumstats g.sumstats --weights weights.pos --weights_dir wgt/ --ref_ld_chr 1KG/EUR. --chr 22 --out chr22.dat`
- CLI (S-PrediXcan, Python): `SPrediXcan.py --model_db_path gtex_v8.db --covariance gtex_v8.cov --gwas_file g.txt --output_file out.csv`
- CLI (S-MultiXcan joint): `SMulTiXcan.py --models_folder mashr_models/ --gwas_folder gwas/ --metaxcan_folder spredixcan_per_tissue/ --output joint.csv`
- CLI (UTMOST cross-tissue): joint test across tissues via UTMOST's per-tissue GBJ / GBJ2 step
- CLI (FOCUS fine-mapping): `focus finemap gwas.sumstats 1KG_EUR focus.db --chr 22 --p-threshold 5e-8 --out chr22.focus`
- CLI (MA-FOCUS multi-ancestry): `focus finemap` with colon-separated per-ancestry sumstats / LD / weights and hyphen-joined ancestry codes in `--locations` (e.g. `38:EUR-EAS-AFR`)

TWAS, cis-eQTL MR, and coloc operate on overlapping evidence: TWAS asks "is the gene's predicted expression associated with the trait?"; cis-eQTL MR asks "does the eQTL effect on expression mediate the trait effect under IV assumptions?"; coloc asks "do the GWAS and eQTL share a causal variant?". Strong causal claims require triangulation, not single-method significance.

## Reference Files

Read the file when the request matches; SKILL.md alone covers model choice, thresholds, triangulation, Common Errors and install.

| File | Read when |
|------|-----------|
| `references/algorithmic-taxonomy.md` | Comparing FUSION, PrediXcan, S-PrediXcan, S-MultiXcan, UTMOST, MOSTWAS, kTWAS, EpiXcan, TIGAR-V2, FOCUS, MA-FOCUS or JEPEG: model, input, strength, failure condition |
| `references/tissue-and-model-selection.md` | Picking the primary tissue (LDSC-SEG / CELLEX / MAGMA), S-MultiXcan vs UTMOST, or cell-type-resolved (sc-eQTL) TWAS |
| `references/failure-modes.md` | Diagnosing LD-induced false positives, tissue mis-specification, ancestry mismatch (with panel table), low-N tissue weights, HLA-region hits, or co-regulated genes |
| `references/fusion-pipeline.md` | Running FUSION per chromosome and the conditional/joint `post_process` step |
| `references/spredixcan-smultixcan.md` | Running per-tissue S-PrediXcan and the S-MultiXcan joint test |
| `references/focus-finemapping.md` | Running FOCUS or MA-FOCUS, building a custom FOCUS DB, `--locations`, Windows paths, the `--prior-prob` scan |
| `references/escalation-and-reconciliation.md` | The standard pipeline misses an expected hit (MOSTWAS, EpiXcan, TIGAR-V2, kTWAS), or FUSION / S-PrediXcan / FOCUS / coloc / MR disagree |
| `references/reporting-and-review.md` | Writing methods for publication or answering reviewer pushback |
| `references/citations.md` | Citing a method or panel |

Scripts (`scripts/`): `fusion_twas.sh` (FUSION per-chromosome TWAS and conditional analysis) and `build_focus_db.py` (FOCUS DB from a weight table). The S-PrediXcan + S-MultiXcan pipeline is `examples/s_predixcan_pipeline.sh`; FOCUS is `examples/focus_finemap.sh`.

## Algorithmic Taxonomy

Tool-by-tool comparison table (model, input, output, strength, fails when): `references/algorithmic-taxonomy.md`. S-PrediXcan and FUSION are mathematically near-identical, so choose by panel availability and ancestry match, not algorithm; the practical differences are the weight panel, LD reference and per-gene heritability threshold.

### PredictDB Model Choice and GTEx Versioning

| Release | Cohort | Status (2026) | When to use |
|---------|--------|----------------|-------------|
| GTEx v8 (838 donors, 49 tissues, 2020) | EUR-dominant (~85%) | Current PredictDB standard | Default; pre-trained MASHR + elastic-net DBs available |
| GTEx v9 (2023) | Expanded harmonization | Not migrated into PredictDB | Do not use until PredictDB rebuilds |
| GTEx v10 (2024 AnVIL release) | Re-aligned to GRCh38 v44 | Limited harmonization; not PredictDB-default | Wait for community-validated weight panels |

**Operational rule:** Use GTEx v8 unless there is an explicit biological reason to deviate (tissue not in v8, ancestry-specific panel preferred). In methods, pin exactly: "GTEx v8 MASHR-EUR, PredictDB release 2022-01".

### MASHR vs Elastic-Net Models

PredictDB ships two cross-validated model families per tissue (Barbeira 2021 Genome Biol 22:49):

| Model | Construction | Per-gene SNP count | When to use |
|-------|--------------|---------------------|-------------|
| MASHR | Cross-tissue posterior mean from DAP-G fine-mapped SNPs | ~10x sparser | Primary discovery; higher per-gene R^2 in most genes; standard for S-MultiXcan |
| Elastic-net | Per-tissue lasso/ridge mix (alpha = 0.5) | Denser | Tissues where MASHR's cross-tissue prior is mis-specified (ovary, testis, isolated-organ traits) |

**Operational rule:** Never mix MASHR and elastic-net within a single S-MultiXcan run; the inter-tissue covariance and condition number assumptions break. Pick one family and apply consistently across all tissues.

## Decision Tree by Experimental Scenario

| Scenario | Recommended workflow | Why |
|----------|---------------------|-----|
| GWAS summary stats only, EUR, single hypothesis tissue (e.g. liver for LDL) | S-PrediXcan with GTEx v8 MASHR-EUR weights, or FUSION with GTEx liver | Standard pre-trained pipeline; minimal compute (`references/fusion-pipeline.md`, `references/spredixcan-smultixcan.md`) |
| GWAS summary stats only, tissue unknown a priori | S-MultiXcan (standard GTEx v8) OR UTMOST (custom panel) -- see S-MultiXcan vs UTMOST table in `references/tissue-and-model-selection.md` | Joint multi-tissue inflates power; tissue prioritisation requires LDSC-SEG separately (`references/tissue-and-model-selection.md`) |
| Multiple TWAS hits at one locus (gene-dense region) | Run TWAS then FOCUS for probabilistic fine-mapping | LD ties co-regulated genes; FOCUS PIP distinguishes likely causal gene (`references/focus-finemapping.md`, `references/failure-modes.md`) |
| Multi-ancestry GWAS (EUR + EAS + AFR) | Per-ancestry S-PrediXcan with matched weights, then MA-FOCUS to combine | Single-ancestry weights miscalibrated in other ancestries; joint fine-mapping shrinks credible set (`references/focus-finemapping.md`) |
| Individual-level genotypes available (UKB) | PrediXcan (full regression) | Allows covariates, interactions, binary outcomes natively |
| Low-N tissue (GTEx N < 100) | Substitute eQTLGen (whole blood, N ~ 31k) OR skip the tissue | Per-gene CV R^2 unstable below N ~ 100; weights overfit (`references/failure-modes.md`) |
| Drug-target prioritisation (TWAS as causal-gene evidence) | TWAS + cis-eQTL MR + coloc + FOCUS triangulation | TWAS alone is associational; triangulation strengthens causal claim |
| Trans-acting / mediator-aware analysis | MOSTWAS | Adds distal trans-mediating SNPs to cis-only models (`references/escalation-and-reconciliation.md`) |
| Rare-variant or non-linear gene effects | kTWAS or TIGAR-V2 | Kernel / non-parametric flexibility (`references/escalation-and-reconciliation.md`) |
| HLA region (chr6:25-35 Mb hg38) | Exclude or use HLA-specific tools | Long-range LD breaks every gene-by-gene method; standard TWAS PIPs not interpretable (`references/failure-modes.md`, HLA region) |
| Splicing-mediated trait (e.g. neuropsych for sQTL) | sTWAS (sQTL-weighted TWAS) using GTEx splicing models | Splicing mediates many GWAS effects; cis-sQTL panels available in PredictDB |
| Cell-type-specific trait | sc-eQTL-based TWAS (e.g. OneK1K, Yazar 2022 Science 376:eabf3041) | Bulk-tissue TWAS averages over cell types; single-cell eQTL recovers cell-type specificity (`references/tissue-and-model-selection.md`) |

## TWAS - MR - Coloc Triangulation

A TWAS-significant gene is associational, not causal. The strongest defensible claim that a gene mediates a GWAS effect comes from triangulating three orthogonal lines of evidence:

| Line | What it tests | Threshold |
|------|---------------|-----------|
| TWAS | Predicted-expression association with trait | S-MultiXcan joint p < 2.3e-6, OR S-PrediXcan per-tissue p < 4.6e-8 (49-tissue Bonferroni), OR per-tissue FDR < 0.05 |
| cis-eQTL MR | Causal effect of expression on trait under IV assumptions (cross-reference causal-genomics/mendelian-randomization) | Wald-ratio or IVW p < 0.05/n_genes; instrument F > 10 |
| Colocalization | Shared causal variant between GWAS and eQTL (cross-reference causal-genomics/colocalization-analysis) | coloc.abf or coloc.susie PP.H4 >= 0.7 |
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
| `DataFrame.pivot() takes 1 positional argument but 4 were given` or `'DataFrame' object has no attribute 'append'` (`focus finemap`, at "Calculating PIPs") | pyfocus 0.802 `finemap.py:1012` / `:188` use APIs removed in pandas 2.0 | Apply the two `finemap.py` patches in Tool Install Notes |
| `focus import` finishes with an empty DB / 0 genes imported | `mygene` or `rpy2` missing (only an ERROR log line); rpy2 also needs R built as a shared library | Install both, or build the DB directly (`scripts/build_focus_db.py`) |
| MA-FOCUS H0 probability dominates | Cross-ancestry heterogeneity at the locus | Run per-ancestry FOCUS separately; do not force joint |
| S-PrediXcan output has effect sizes much larger than expected | `sdY` proxy mis-specified; standardised vs unstandardised mismatch | Confirm GWAS Z and beta scale; rerun with `--additional_output` |
| `InvalidArguments: Specify either cutoff_ratio or cutoff_threshold` (S-MultiXcan) | `--cutoff_condition_number` (or `--cutoff_ratio`/`--cutoff_threshold`) omitted; not optional in installed versions | Always pass `--cutoff_condition_number 30` |
| A comma-delimited filter over S-MultiXcan's `--output *.csv` silently returns 0 rows (no error) | `SMulTiXcan.py` hard-codes tab-separated output regardless of the output filename's extension | Filter with `awk -F'\t'` (see `examples/s_predixcan_pipeline.sh`), not `-F','` |
| `TypeError: read_csv() got an unexpected keyword argument 'delim_whitespace'` (`focus`/`focus finemap`) | Unpinned pandas>=2.2 | Install pins and patch in Tool Install Notes |
| `AttributeError: module 'numpy' has no attribute 'warnings'` (`focus`/`focus finemap`) | Removed `np.warnings` alias; pinning pandas does not help | Patch after install (Tool Install Notes) |
| `Please specify independent regions location or default regions with '37:EUR', etc.` (`focus finemap`) | `--locations` omitted | Always pass it (see `references/focus-finemapping.md`) |
| `focus finemap` reports an implausible population count (e.g. "Detecting 2 populations" from one file) on Windows | Absolute Windows paths get colon-split by pyfocus's multi-ancestry parser | Use relative paths (see `references/focus-finemapping.md`) |
| `FileNotFoundError: ...bim` on the `ref` argument to `focus finemap` | LD reference path templated like FUSION's per-chromosome `--ref_ld_chr` (e.g. `.../chr` with no matching file); pyfocus passes it straight to `pandas_plink.read_plink()`, which needs an exact prefix or a real glob, not a chr-substitution convention | Point at one PLINK fileset covering every chromosome analyzed, or a `pandas_plink` glob (`prefix/chr*.bed`); `--chr`/`--locations` filter after loading, they do not select the file (see `references/focus-finemapping.md`) |
| A PIP filter against a `pip` column silently returns nothing (or the whole row) | Installed pyfocus 0.802 never writes a plain `pip` column -- output is `pips_pop1` (population-indexed even for a single population) or `pips_me` for the cross-ancestry marginal PIP | Filter on `pips_pop1` (single-ancestry) or `pips_me` (MA-FOCUS); see `examples/focus_finemap.sh` |
| `Error in wgt.matrix[qc$flip, ] : incorrect number of dimensions` or `non-conformable arguments` (`FUSION.post_process.R`) | Single-SNP ("top1") weight models: `wgt.matrix[m.keep,]` and `genos$bed[,m[m.keep]]` (lines ~168/170 and ~251/254) drop a 1-row/1-column matrix to a bare vector without `drop=FALSE` | Patch your local `fusion_twas` clone: add `,drop=FALSE` to both subsetting lines at both locations, or filter single-SNP genes out of `--input` before running `post_process.R` |

## Tool Install Notes

- **FUSION**: gusevlab.org/projects/fusion or `git clone https://github.com/gusevlab/fusion_twas`. R scripts; needs plink (PLINK 1.9) and Rscript; `install.packages(c('plink2R', 'optparse', 'glmnet', 'methods', 'RColorBrewer', 'GBJ'))` (GBJ is for the omnibus test). Pre-computed weights for GTEx v7/v8, CMC, YFS, METSIM, NTR, MESA at the same site.
- **MetaXcan / S-PrediXcan / S-MultiXcan**: `git clone https://github.com/hakyimlab/MetaXcan` (not on PyPI). Python scripts in `software/`. Compatible with Python 3.9-3.11.
- **PredictDB models**: predictdb.org. GTEx v8 elastic-net (single-tissue) and MASHR (cross-tissue posterior) databases for EUR; multi-ethnic panels emerging.
- **UTMOST**: `git clone https://github.com/Joker-Jerome/UTMOST`. Python. Needs precomputed cross-tissue weights or training pipeline.
- **FOCUS**: `pip install pyfocus "pandas<2.2" "setuptools<81"` -- checked against pyfocus 0.802 (pandas 2.1.4, numpy 1.26.4, setuptools 80.10.2), 2026-09-21. A bare `pip install pyfocus` is non-functional and the pins alone are not enough; four upstream problems, the first two fixed by the pins and the last two by the patch below (unfixed upstream as of 2026-09-21):
  - `delim_whitespace`: pandas>=2.2 removed it (used in `gwas.py`/`exprref.py`/`ldref.py`/`convert.py`), hence `pandas<2.2`.
  - `pkg_resources`: imported by pyfocus but not declared as a dependency, and gone from setuptools>=81, hence `setuptools<81`.
  - `np.warnings`: removed from every numpy that has Python 3.11/3.12 wheels (`ldref.py`, `exprref.py`).
  - `DataFrame.pivot(...)` with positional args (`finemap.py:1012`) and `DataFrame.append` (`finemap.py:188`): both removed in pandas 2.0, and both crash at "Calculating PIPs", after every earlier stage has passed.

  ```bash
  PYFOCUS_DIR=$(python -c "import pyfocus, os; print(os.path.dirname(pyfocus.__file__))")
  sed -i "1i import warnings" "$PYFOCUS_DIR/data/exprref.py"
  sed -i "s/np\.warnings\.catch_warnings/warnings.catch_warnings/;s/np\.warnings\.filterwarnings/warnings.filterwarnings/" \
      "$PYFOCUS_DIR/data/ldref.py" "$PYFOCUS_DIR/data/exprref.py"
  sed -i 's/\.pivot("model_id", "attr_name", "value")/.pivot(index="model_id", columns="attr_name", values="value")/' "$PYFOCUS_DIR/finemap.py"
  sed -i 's/^\( *\)df = df\.append(null_dict, ignore_index=True)/\1df = pd.concat([df, pd.DataFrame([null_dict])], ignore_index=True)/' "$PYFOCUS_DIR/finemap.py"
  ```

  Verified 2026-09-21 from a fresh venv with this pin + patch: `focus finemap` against a real pyfocus-schema DB (built as under "FOCUS Probabilistic Fine-Mapping") runs to completion and returns `pips_pop1 = 1` for the planted true gene and 3.4e-08 for `NULL.MODEL`. CLI `focus`. Pre-built DBs for GTEx v7/v8 panels at github.com/bogdanlab/focus. `focus import` additionally needs `mygene` + `rpy2` -- see the FOCUS section.
- **MA-FOCUS**: `git clone https://github.com/mancusolab/ma-focus && cd ma-focus && pip install .` (no PyPI release; install from source). Same CLI as single-ancestry FOCUS -- `focus finemap` -- with colon-separated per-ancestry sumstats / LD / weight DBs and paired ancestry codes in `--locations`.
- **TIGAR-V2**: `git clone https://github.com/yanglab-emory/TIGAR`. Python + R hybrid; ships with example data.
- **MOSTWAS**: `git clone https://github.com/bhattacharya-a-bt/MOSTWAS`. R package; install via `devtools::install_github`.
- **EpiXcan**: Bitbucket roussoslab/epixcan. Workflow-style pipeline.

## Related Skills

- causal-genomics/fine-mapping - Variant-level credible sets feeding FOCUS gene-level fine-mapping
- causal-genomics/colocalization-analysis - Coloc PP.H4 triangulation with TWAS hits
- causal-genomics/mendelian-randomization - cis-eQTL MR triangulation; drug-target prioritisation
- causal-genomics/effector-gene-prioritization - Downstream gene mapping from TWAS candidate sets
- causal-genomics/proteome-mr-drug-target - Drug-target triangulation using protein-level evidence
- causal-genomics/pleiotropy-detection - Distinguishing horizontal pleiotropy from mediated TWAS signal
- causal-genomics/mediation-analysis - Downstream gene-mediated trait effects given TWAS hits
- population-genetics/association-testing - Upstream GWAS summary statistic generation
- population-genetics/linkage-disequilibrium - LD reference panel construction for FUSION / FOCUS
- differential-expression/deseq2-basics - Generating eQTL count data for custom prediction-weight training
- single-cell/preprocessing - Cell-type-resolved eQTL panels for sc-TWAS
- workflows/gwas-pipeline - End-to-end GWAS pipeline producing TWAS input
- variant-calling/variant-annotation - Functional annotation of TWAS / FOCUS top variants
