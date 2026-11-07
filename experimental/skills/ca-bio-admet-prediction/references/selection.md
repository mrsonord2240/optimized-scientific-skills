## ADMET Model Taxonomy

| Tool | Endpoints | Architecture | Uncertainty | Access | Fails when |
|------|-----------|--------------|-------------|--------|------------|
| ADMETlab 3.0 | 119 reported features: 77 prediction models, 34 computed properties, 8 rules | Multi-task DMPNN + descriptors for modeled endpoints | Evidential uncertainty for modeled endpoints | Web service; hosted API documented by the authors | Outside training distribution; metals; macrocycles |
| ADMET-AI | TDC-derived ADMET tasks; inspect installed model metadata | Chemprop D-MPNN | Inspect version-specific outputs; do not assume calibrated uncertainty | Python package | v2 package predictions differ from the v1 paper/server |
| DeepChem MolNet | Dataset-dependent tasks including Tox21, ToxCast, and ClinTox | Model-dependent | Model-dependent | Python package | Coverage and uncertainty depend on the selected dataset/model |
| pkCSM | Service-defined ADMET endpoints | Graph signatures + ML | Inspect current service output | Web service | Applicability domain and service contract must be checked |
| SwissADME | Physchem, pharmacokinetics, drug-likeness, and medchem outputs | Published models and rules | None advertised | Web service (no public API) | Automated access is restricted by its terms |
| ProTox-3.0 | 61 toxicity models/endpoints | RF/DNN + fingerprints, similarity, and pharmacophore methods | Confidence score | Web service / sample API | Toxicity only; reports LD50 and toxicity class |
| ADMETpredictor (Simulations Plus) | ~140 | Proprietary | Per-prediction | Commercial | License cost |
| FAF-Drugs4 | filters | Rule-based | None | Web | Static rules |
| chemprop (in-house) | User-defined | D-MPNN ± descriptors | Ensemble and other estimators; optional calibration | Python package | Requires suitable training and calibration data |

**Decision:** For a reproducible offline first pass, use the bundled **ADMET-AI 2.x** route with the input and applicability-domain gate below. Use **ADMETlab 3.0** only after obtaining its current official API contract when its broader hosted coverage and modeled-endpoint uncertainty are needed. For a sufficiently large, relevant in-house endpoint dataset, benchmark a **chemprop D-MPNN**, descriptors, and simpler baselines under a deployment-relevant split rather than assuming a universal sample-size threshold. Shan et al. (2022) reported an AUC of 0.956 for a D-MPNN combined with 206 MOE descriptors on their random-split hERG benchmark.

## Decision Tree by Scenario

| Scenario | Workflow | Reasoning |
|----------|----------|-----------|
| Library triage, no in-house data | Bundled ADMET-AI 2.x batch + reference-similarity gate | Reproducible offline route; use its point predictions only after the gate |
| Need ADMETlab's 119 reported platform features | Current official ADMETlab 3.0 API workflow + `load_admetlab_results` | Hosted coverage and modeled-endpoint uncertainty; retain task ID and uncertainty columns |
| Single endpoint, adequate in-house data | Benchmark chemprop D-MPNN, descriptors, and simpler baselines | Select by prospective or deployment-relevant validation |
| Need calibrated probabilities | chemprop with ensemble + Platt | Native deep learning rarely calibrated |
| FDA / regulatory submission | OECD-compliant QSAR with AD | See OECD principles below |
| Quick annotation for VS | Lipinski, Veber, and QED reported separately | Rank or annotate; do not impose a universal QED gate |
| BBB penetration | Simple screen: TPSA <= 90, MW <= 500, HBD <= 3 | Repository heuristic; not the six-factor CNS MPO |
| Cardiotox liability | ADMETlab hERG + ProTox-3.0 cardiotoxicity + literature check | Both hosted endpoints model hERG blockade; compare applicability domains and assay definitions |
| Drug-drug interaction (CYP) | CYP1A2/2C9/2C19/2D6/3A4 inhibitor + substrate | Standard set of 5 CYPs |

## OECD QSAR Principles (5 Pillars)

For regulatory-grade ADMET QSAR (REACH, ECHA, FDA submissions), models must satisfy:

1. **Defined endpoint** -- specific bioassay, units, conditions
2. **Unambiguous algorithm** -- reproducible model + code
3. **Defined applicability domain (AD)** -- where the model is valid
4. **Appropriate statistical validation** -- external test set, cross-validation
5. **Mechanistic interpretation** -- biological / chemical rationale

For non-regulatory work, AD assessment is still critical. The OECD's *applicability domain* is the workhorse: predictions outside the AD are unreliable, but operational AD measures (leverage, kNN, conformal prediction) often disagree.

## Applicability Domain Methods

| Method | Definition | Flags out-of-AD when |
|--------|-----------|----------------------|
| kNN distance | Mean distance to k nearest neighbors in training set | > training-set distribution P95 |
| Leverage (Williams) | Hat-matrix diagonal | > 3p/n (p = features, n = compounds) |
| Density (KDE on PCA) | Density in feature space | < density of training set P5 |
| Conformal prediction | Per-prediction confidence interval | Interval > tolerance |
| Bayesian variance | Ensemble or MC-dropout variance | > training-set variance P95 |

For deep-learning ADMET, **conformal prediction** can provide calibrated prediction sets or intervals when its exchangeability and calibration assumptions are appropriate (McShane et al. 2024).
