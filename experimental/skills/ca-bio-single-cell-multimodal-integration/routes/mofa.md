# MOFA+ factors

```bash
python scripts/mofa_muon.py multiome.h5mu mofa_out --n-factors 15
```

Input: a `.h5mu` with each modality already in the representation its likelihood needs: scaled RNA gaussian, binarized ATAC bernoulli, counts poisson. The likelihoods must match the data. Output: `mofa_out.h5mu` with `obsm['X_mofa']` and `mofa_out.hdf5` (model, loadings, variance explained per modality).

- MOFA+ is for interpretation: shared versus modality-specific axes of variation. Do not use it for clustering or denoising, and it does not correct batch within a view.
- Interpret factors through loadings and per-modality variance explained.

Done when `mofa_out.hdf5` exists and the answer reports variance explained per factor and modality.
