## Predictor Taxonomy

| Family | Architecture | Output | Fails when |
|--------|--------------|--------|------------|
| Context-aware CNN | 10 kb dilated ResNet | Per-position donor/acceptor probability | Long-range (>5 kb) regulatory effects; tissue-specific events |
| Tissue-aware CNN/transformer | Same arch + multi-tissue training | Per-tissue ΔPSI | Tissue not in training set; novel cell types |
| Modular per-region CNN | Separate sub-models for 5'ss/3'ss/exon/intron | Calibrated quantitative ΔPSI | Atypical events; complex multi-junction effects |
| Foundation transformer | Pretrained on broad genomic context | Splice probability or ΔPSI | New tools; less battle-tested |
| Empirical lookup | Public RNA-seq event database | Top-N most likely mis-splicing outcomes | Variant types not represented in training cohorts |
| Composite score | Blend of multiple predictors | Single scaled score | When component predictors disagree internally |

## Tool Selection Matrix

| Tool | Best for | Output | When to use | Fails when |
|------|----------|--------|-------------|------------|
| SpliceAI | Clinical screening; canonical splice site disruption | Delta score 0-1 | Default for ACMG variant classification | Tissue-specific events; effect site farther than -D from the variant |
| Pangolin | Tissue-aware predictions | Max gain/loss over 4 tissue models (per-tissue via Python) | When disease tissue is known (brain, heart, liver, testis) | Tissue not in 4-tissue training set; `-m True` (see Pangolin) |
| MMSplice | Quantitative ΔPSI | Δlogit_psi | Research where calibrated effect-size matters | Atypical events outside cassette-exon model; variants far from an exon |
| SpliceTransformer | 2024+ benchmark improvements | Splice score + 15 tissue-specificity flags (SNVs) | When transformer foundation models outperform CNN on benchmark variant sets | New (2024); limited clinical adoption |
| CI-SpliceAI | Second SpliceAI-family opinion | SpliceAI-style delta scores | Extended windows; cross-check | Did not flag GLA c.639+919G>A here (0.05) |
| SpliceVault | Empirical mis-splicing outcome | Top-N events at the affected splice site | Predicting consequence (skip vs cryptic) of canonical-disrupting variants | Variants not at catalogued splice sites |
| CADD-Splice | Single composite score | Scaled C-score (PHRED) | Clinical pipelines wanting one number | When knowing which sub-component drove the score is needed |

Methodology evolves; verify benchmarks (Smith & Kitzman 2023 *Genome Biol* 24:294; You et al 2024 *Nat Commun*) and ClinGen SVI splicing recommendations before reporting any interpretation. Concordance is a prioritization signal, not proof: even three concordant predictors gave a false-positive splice-disruption signal for ClinVar likely-benign CFTR c.2909-10T>G in the checked panel. Check ClinVar/gnomAD, disease mechanism and RNA evidence; discordance also flags RNA-validation candidates.

## Decision Tree by Use Case

| Use case | Recommended approach |
|----------|----------------------|
| Single variant, ACMG-style evidence | SpliceAI default 50nt + ClinGen SVI 2023 thresholds |
| Tissue-specific question (brain disease, cardiomyopathy) | SpliceAI + Pangolin per-tissue (`pangolin_tissue.py`) or SpliceTransformer |
| Unsolved Mendelian case (suspect deep-intronic) | SpliceAI `-D 500` on every intronic candidate, both pseudoexon ends reported; SpliceVault for the consequence |
| VUS panel screening | SpliceAI + Pangolin + MMSplice concordance (`splice_parsers.py`) |
| Predict consequence of canonical-disrupting variant | SpliceVault top-N empirical events |
| Branchpoint variant suspected | Predictors are weak here; RNA validation (see Branchpoint Variant Detection) |
| Splice-switching ASO design | Checklist only (see ASO section); no executable here |
| Validate predicted splice change in patient | RNA-seq + FRASER2 (see outlier-splicing-detection) |
| Pseudoexon prediction in deep intron | SpliceAI extended window; CI-SpliceAI as second opinion; require RNA validation |

## ClinGen SVI 2023 Framework

The ClinGen Sequence Variant Interpretation (SVI) splicing subgroup (Walker 2023 *Am J Hum Genet*) extended the ACMG/AMP 2015 framework with explicit splice-prediction rules.

| Evidence code | Threshold | Notes |
|----------------|-----------|-------|
| **PP3** (supporting pathogenic) | SpliceAI delta >= 0.20 | ClinGen SVI: apply at **supporting** weight (not standalone) |
| **BP4** (supporting benign) | SpliceAI delta <= 0.10 | ClinGen SVI: apply at **supporting** weight |
| **PVS1** (very strong null) | Canonical +/-1, +/-2 site disruption with predicted LoF + NMD | Requires gene where LoF is established mechanism (Abou Tayoun 2018 *Hum Mutat* PVS1 decision tree) |
| **PS3 / BS3** (functional) | RNA evidence (RT-PCR, RNA-seq, minigene) | Supersedes computational evidence |
| *(no code)* | delta 0.10-0.20, or **no score** (`.`, no `SpliceAI=` tag) | Inconclusive / `not_scored`: a missing score is never BP4 |

Boundaries are inclusive as published: 0.10 -> BP4, 0.20 -> PP3. SpliceAI scores have two decimals, so exact 0.20/0.50/0.80 values are common.

**Operational rules:** Computational evidence (PP3/BP4) is *supporting*, not standalone. Higher SpliceAI cutoffs (0.5, 0.8) increase precision but are the tool's own tiers (Jaganathan 2019), NOT ClinGen-endorsed evidence-strength upgrades — reaching moderate/strong requires functional/RNA evidence (PS3/BS3), not a higher SpliceAI score alone. Splicing variants benefit from concordance across SpliceAI + Pangolin + MMSplice. RNA validation supersedes prediction. Always log SpliceAI version, distance window, and reference transcript. SpliceAI alone is **not sufficient** for PVS1 (canonical site disruption also needs gene-level LoF context) or for non-canonical positions.
