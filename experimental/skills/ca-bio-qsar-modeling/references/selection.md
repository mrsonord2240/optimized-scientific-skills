## Model Taxonomy

| Model | Architecture | Use case | Fails when |
|-------|--------------|----------|------------|
| Random Forest + ECFP4 | Classical baseline | Small-data comparison, interpretability | May miss signal not represented by the fingerprint |
| chemprop D-MPNN | Directed message passing | Graph-learning candidate to benchmark | Can overfit when data are sparse or biased |
| chemprop D-MPNN + RDKit 2D | Hybrid graph + descriptors | Useful hybrid baseline; compare on the same split | Diminishing returns at large data |
| MolFormer | SMILES transformer | Large public training data benefit | Compute overhead; OOD risk |
| Uni-Mol | 3D-aware transformer | 3D-relevant endpoints (binding) | Requires 3D conformers |
| ChemBERTa-2 | SMILES transformer pretrained on up to 77M molecules | SMILES language-model baseline | Fine-tuning benefit is endpoint- and split-dependent |
| Gaussian Process + ECFP4 | Probabilistic | Active learning; uncertainty | O(N^3) scaling |
| MultiTask DNN | Joint training | Multiple endpoints | Data must overlap |

**Decision:** Compare a fingerprint-based baseline with chemprop under the same split and endpoint. Add a pretrained transformer or 3D model only when its representation, compute cost, and validation design fit the deployment question; dataset size alone does not determine the winner.

## Decision Tree by Scenario

| Dataset context | Endpoint type | Model to benchmark |
|--------------|---------------|-------|
| Sparse labels or few independent series | Regression / classification | Regularized fingerprint baseline; quantify instability and avoid unsupported deployment |
| Multiple scaffold groups with adequate labels | Regression / classification | Fingerprint baseline plus chemprop on identical splits |
| Large public or internal training collection | Regression / classification | Benchmark chemprop and a relevant pretrained representation |
| Multi-task | Related endpoints (CYP3A4, CYP2D6, etc.) | chemprop MultiTask |
| 3D-relevant | Binding, conformer-dependent | Uni-Mol with conformer ensemble |

## OECD 5 Principles

The OECD principles were agreed in 2004; the 2007 guidance explains their application:

1. **Defined endpoint**: specific bioassay, units, threshold definitions
2. **Unambiguous algorithm**: reproducible code, fixed random seeds, version-pinned dependencies
3. **Defined applicability domain (AD)**: where the model is valid
4. **Appropriate measures of goodness-of-fit, robustness, and predictivity**: external test set and suitable validation
5. **Mechanistic interpretation, if possible**: biological/chemical rationale where available

For non-regulatory QSAR, all 5 still good practice; especially **AD definition** is critical.

## Applicability Domain Methods

| Method | Definition | Pro | Con |
|--------|-----------|-----|-----|
| **Ensemble variance** | Variance across N-model ensemble predictions | Supported by `chemprop predict --uncertainty-method ensemble` when multiple model paths are supplied | Assumes useful ensemble diversity; not calibrated coverage; variance is in squared endpoint units |
| kNN distance | Mean Tanimoto to k nearest in training | Easy to interpret | Doesn't account for label distribution |
| Leverage | Hat-matrix diagonal | Statistical on low-dimensional, well-conditioned continuous descriptors | Do not apply directly to sparse ECFP bits; reduce/select descriptor space first and reject a calculation with any `h < 0` or `h > 1` |
| KDE on PCA | Density in feature space | Captures multivariate structure | Density choice subjective |
| Mahalanobis distance | Covariance-aware distance | Theoretically motivated in a low-dimensional, well-conditioned descriptor/PCA space | Reduce/select descriptor space first; unstable on high-dimensional fingerprints |
| Conformal prediction | Per-prediction interval or set | Finite-sample marginal coverage under exchangeability | Requires a calibration design and compatible predictor; evaluate held-out empirical coverage after scaffold/time splits |
| Bayesian / MC-dropout | Posterior or dropout variance | Direct uncertainty | Computational cost |
| Tanimoto coverage | At least 1 NN within threshold | Practical | Threshold subjective |

Ensemble disagreement is one useful uncertainty diagnostic, not a formally defined applicability domain or calibrated coverage guarantee. If using a threshold such as a training-distribution percentile, label it as a project-defined heuristic and validate it prospectively.
