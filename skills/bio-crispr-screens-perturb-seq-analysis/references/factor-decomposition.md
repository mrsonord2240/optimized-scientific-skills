## Factor-Based Analysis

For complex perturbation responses, decompose the per-cell perturbation effect into shared latent factors:

```bash
# FR-Perturb ("Factorize-Recover") decomposes perturbation effects into shared factors.
# It is NOT part of pertpy: it is a standalone CLI from douglasyao/FR-Perturb
# (Yao et al. 2023 Nat Biotechnol). Install is clone + conda env from its environment.yml
# (Python 3.8; its `python-spams` dependency is only on the `defaults`/`anaconda` channel), so
# run it outside your analysis env. Flags below checked against run_FR_Perturb.py's argparse
# (2026-09-22) -- there is no `--input`/`--perturbations` form.
# Perturbation status either as a cell x perturbation matrix:
python run_FR_Perturb.py --input-h5ad counts.h5ad \
    --input-perturbation-matrix perturbations.txt --out results/prefix
# ...or from the .obs column (use --perturbation-delimiter for multi-perturbation cells):
python run_FR_Perturb.py --input-h5ad counts.h5ad \
    --perturbation-column-name <obs column> --out results/prefix
# Useful additions: --control-perturbation-name, --covariates, --output-factor-matrices.
```
