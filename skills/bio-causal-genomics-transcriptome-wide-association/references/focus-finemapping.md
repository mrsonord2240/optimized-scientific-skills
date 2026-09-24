# FOCUS Probabilistic Fine-Mapping

**Goal:** Resolve a locus with multiple co-significant TWAS genes into a probabilistic causal-gene credible set.

**Approach:** Build (or download) a FOCUS gene-prediction database matching the TWAS weight panel; run `focus finemap` with the TWAS sumstats and ancestry-matched LD reference; report PIPs and credible sets.

```bash
# Install: see SKILL.md Tool Install Notes -- the pip pins plus a post-install patch are required to run at all.
# Pre-built FOCUS DBs for FUSION/PrediXcan weights live at github.com/bogdanlab/focus

# Convert FUSION TWAS output to FOCUS sumstat format if needed
# FOCUS expects: CHR SNP BP A1 A2 Z P from the underlying GWAS sumstats (NOT TWAS Z)

focus finemap \
    gwas.sumstats \
    1000G_EUR_all \
    focus_gtex_v8_whole_blood.db \
    --p-threshold 5e-8 \
    --tissue Whole_Blood \
    --locations 38:EUR \
    --out gwas_focus_whole_blood
# Output: per-gene PIP, credible-set membership flag, and locus-level group probability
```

**LD reference path is a single PLINK fileset, not per-chromosome like FUSION.** `focus finemap` passes the `ref` argument straight to `pandas_plink.read_plink()`: give it one bfile prefix (`1000G_EUR_all.bed/.bim/.fam`) covering every chromosome you will analyze, or a real glob (`1000G_EUR/chr*.bed`) if the panel is split into per-chromosome files. `--chr`/`--locations` only filter the SNPs already loaded -- they do not select which file to open. A chr-templated prefix with no matching file (e.g. `1000G_EUR_chr`, mirroring FUSION's `--ref_ld_chr`) 404s with `FileNotFoundError: ...bim` (confirmed 2026-09-21).

`--locations` is required on every `focus finemap` call in installed pyfocus 0.802, including single-ancestry runs -- omitting it crashes with `Please specify independent regions location or default regions with '37:EUR', etc.` (confirmed 2026-09-19; there is no "uses default LD blocks" fallback). Use `38:EUR` for GRCh38-aligned panels (e.g. GTEx v8 PredictDB/MASHR) or `37:EUR` for GRCh37-aligned panels; non-EUR/multi-ancestry codes follow the same `build:ANC1-ANC2-...` syntax used by MA-FOCUS below.

**Windows path caveat:** pyfocus 0.802 splits the `gwas`/`ref`/`weights` positional arguments on `:` to detect multi-ancestry input (see MA-FOCUS below) -- this also splits the drive-letter colon in any absolute Windows path (e.g. `F:/data/gwas.sumstats` is parsed as 2 populations), even for a single-ancestry call. Confirmed 2026-09-19: an absolute `F:/...` path made `focus finemap` mis-detect "2 populations" from one file. Use relative paths (from the working directory) for every `focus finemap` argument on Windows, or run under WSL/Linux.

For a custom prediction-weight panel without a pre-built FOCUS DB, construct one from FUSION weights:

```bash
focus import custom_panel.pos fusion --tissue Whole_Blood --output custom_focus
```

`focus import` needs `mygene` and `rpy2` installed (`pip install mygene rpy2`), and rpy2 needs R built as a shared library. Without them it logs only an ERROR-level message and leaves an empty DB (0 genes imported) instead of raising -- check the log and the DB's gene count. Where rpy2 is unavailable (e.g. the Windows R here), build the DB directly with pyfocus's own schema from a table `panel.tsv` (columns `gene chrom txstart txstop snp pos a1 a0 weight`, one row per gene-SNP weight); `focus finemap` reads the result as a normal DB (checked 2026-09-21, pyfocus 0.802, SQLAlchemy 2.0.54):

```bash
python <skill-dir>/scripts/build_focus_db.py panel.tsv custom_focus.db --ref-name custom_panel --tissue Whole_Blood
```

The script writes one pyfocus model per gene (`method="top1"`); pass each gene's real CV R2 and p as optional `cv_r2` / `cv_r2_pval` columns of `panel.tsv` (defaults 0.1 and 1e-4 are placeholders).

FOCUS PIPs depend on the per-locus prior probability that any gene is causal. Report a sensitivity scan over `--prior-prob`:

| Setting | Use case |
|---------|----------|
| `--prior-prob 1e-3` (default) | Standard genome-wide TWAS; matches Mancuso 2019 |
| `--prior-prob 1e-4` (conservative) | High-prior gene-dense locus where most genes are not causal |
| `--prior-prob 1e-2` (liberal) | Pre-prioritised candidate region where one gene is expected |

Cite the Mancuso 2019 supplement for the sensitivity-scan protocol. MA-FOCUS extends with per-ancestry weights and is invoked via the same `focus finemap` CLI -- multi-ancestry mode is signaled by colon-separated per-ancestry sumstats, LD references, and weight DBs, plus paired ancestry codes in `--locations`:

```bash
focus finemap \
    eur.sumstats.tsv.gz:eas.sumstats.tsv.gz:afr.sumstats.tsv.gz \
    1000G_EUR_all:1000G_EAS_all:1000G_AFR_all \
    focus_eur.db:focus_eas.db:focus_afr.db \
    --chr 22 --locations 38:EUR-EAS-AFR \
    --out gwas_ma_focus
```

MA-FOCUS assumes a shared causal gene across ancestries; gene-specific heterogeneity (e.g. an ancestry-specific eQTL) violates the assumption and produces inflated H0/heterogeneous-group probability.
