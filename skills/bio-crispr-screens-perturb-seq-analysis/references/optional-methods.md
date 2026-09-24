## Optional-method environment boundaries

The default, executable Python path is `assign_sgrna.py` + `mixscape_filter.py` + Pertpy/PyDESeq2. It is the fallback when either optional method below is absent. Do not present an unavailable optional dependency as if it ran.

### FR-Perturb

FR-Perturb is a standalone project, not a Pertpy feature. Its upstream environment is Python 3.8 with `python-spams` from the `defaults`/`anaconda` channel. Keep it isolated from a modern (for example Python 3.12) Pertpy environment:

```bash
git clone https://github.com/douglasyao/FR-Perturb.git
cd FR-Perturb
conda env create -f environment.yml
conda activate fr-perturb
python run_FR_Perturb.py --help
```

If that environment cannot be created, use the supplied Pertpy/PyDESeq2 per-perturbation DE workflow and report that factor decomposition was not run.

### Seurat

Seurat/Mixscape is optional R tooling. It is not required for the supplied Python Mixscape and multiome scripts. Install and verify it only in an R environment intended for Seurat:

```r
install.packages("Seurat")
stopifnot(requireNamespace("Seurat", quietly = TRUE))
packageVersion("Seurat")
```

If Seurat is unavailable, use `scripts/mixscape_filter.py`; for RNA+ATAC use `scripts/multiome_differential.py`. Record the fallback in the analysis provenance.
