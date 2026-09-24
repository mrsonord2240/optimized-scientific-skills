# Hit-Calling Method Catalog

Read when choosing between methods or explaining why a method fits a design. The decision tree in SKILL.md gives the pick; this is the comparison behind it.

## Statistical Models Compared

| Method | Year | Statistical model | Tests | Best for | Fails when |
|--------|------|-------------------|-------|----------|------------|
| MAGeCK RRA | 2014 | NB per-sgRNA -> alpha-RRA per gene | Two-sided | General two-condition | >40% guides change (median norm breaks); time course; cancer-line CN |
| MAGeCK MLE | 2015 | NB GLM with design matrix; per-gene beta | Wald per condition | Multi-condition / time course | Cell-line specific essentiality; CN bias |
| BAGEL2 | 2021 | Bayes factor from log-likelihood ratio | Essential vs non-essential | Essentiality classification | Non-essentiality screens; drug screens |
| drugZ | 2019 | Bidirectional Z-score on guide-level LFC | Sensitizer vs suppressor | Drug-modifier / chemogenomic | Essentiality (no biological prior); time-course |
| JACKS | 2019 | Variational Bayes: LFC = gene * efficacy | Per-gene posterior | Multi-screen joint, library calibration | Single screen; cross-chemistry |
| Chronos | 2021 | Cell-population dynamics ODE + NB | Gene effect adjusted for screen quality | Cancer-line panels, longitudinal | Single screen; non-cancer applications |
| CERES | 2017 | Nonlinear model decoupling CN-bias from gene effect | Per-gene effect | Cancer-line panel with CN profile | Superseded by Chronos at DepMap |

## RRA vs MLE Within MAGeCK

| Property | RRA (`mageck test`) | MLE (`mageck mle`) |
|----------|----------------------|---------------------|
| Conditions supported | 2 | Multiple (design matrix) |
| Statistical test | Robust rank aggregation | Wald on beta from NB GLM |
| Output | neg/pos score, FDR per direction | beta per condition |
| sgRNA efficiency | Not modeled (optional fixed input) | Modeled via `--sgrna-efficiency` |
| Outlier robustness | High (rank-based) | Lower (likelihood-based) |
| Best for | Standard 2-condition screen | Time course, drug screen, multi-cell-line, paired |
| Speed | Fast | Slow (per-gene optimization) |

## Algorithmic Taxonomy: Why Each Was Built

| Method | Designed to solve |
|--------|--------------------|
| MAGeCK RRA | First robust statistical framework for CRISPR-screen ranking; alpha-RRA borrowed from RRA in microarray meta-analysis |
| MAGeCK MLE | Extend MAGeCK to multi-condition; explicit beta scores allow direct LFC interpretation |
| BAGEL2 | Reference-set-anchored Bayesian classification; precision-recall calibrated; tumor-suppressor sensitivity (BAGEL1 was uni-directional) |
| drugZ | Drug-modifier screens have low effect sizes and need bidirectional sensitivity; STARS/MAGeCK miss synthetic-lethal hits |
| JACKS | Sample-size reduction via library-shared efficacy; library calibration as side product |
| Chronos | DepMap-scale (1000+ cell lines, billions of cell-divisions) needs population-dynamics model; CN bias + screen quality first-class |
| CERES | First to formally decouple CN from gene effect at DepMap scale; superseded but historically important |

