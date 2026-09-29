
## Version Compatibility

Runnable examples tested with requests 2.31+ and the public coordinate-based GeneBe endpoint. AutoPVS1 (Xiang 2020) and InterVar 2.2+ are documented external tools, not shipped integrations. ACMG/AMP Bayesian point system is Tavtigian 2018 *Genet Med* / 2020 *Hum Mutat*. Pejaver 2022 *AJHG* and Bergquist 2025 *Genet Med* provide PP3/BP4 calibration. ClinGen Splicing Subgroup 2023 (Walker *AJHG*) governs splicing evidence.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. VCEP-specific CSpec rules override default ACMG application. Use the current CSpec REST service at `https://cspec.genome.network/cspec/srvc`, resolve a gene's available versions via `/cspec/Gene/id/{GENE}/SequenceVariantInterpretation/version`, then record the selected specification version and access date. The former `/cspec/ui/svi/all` GET route is stale.

# ACMG/AMP Framework Reference

**'Classify this variant per ACMG/AMP'** -> Apply 28-criterion framework using Tavtigian point system; gate on ClinGen SVI specifications and VCEP-specific overrides; assign P / LP / VUS / LB / B classification with evidence trail.

- Python (automated): GeneBe API `https://api.genebe.net/cloud/api-public/v1/variant`
- Python (rule-based): InterVar -> `python InterVar.py -i input.vcf -b hg38 --table_annovar table_annovar.pl`
- Web tools: VarSome (commercial), Franklin/Genoox (commercial), ClinGen VCI (gold standard for SVI)
- Citation: Richards 2015 *Genet Med* 17:405 (original framework); Tavtigian 2020 *Hum Mutat* 41:1734 (point system)

## The Tavtigian Bayesian Point System: The Engine Inside All Modern Classifiers

Richards 2015 specified 28 criteria with strength labels (Supporting / Moderate / Strong / Very Strong); combination rules produced P / LP / VUS / LB / B. **Tavtigian 2018/2020 demonstrated this framework is mathematically a Bayesian classifier and proposed the naturally-scaled point system that every modern automated classifier implements:**

| Strength | Points | Odds of pathogenicity |
|----------|--------|----------------------|
| Supporting | 1 | 2.08:1 |
| Moderate | 2 | 4.33:1 |
| Strong | 4 | 18.7:1 |
| Very Strong | 8 | 350:1 |

Benign codes are negative-signed. Final classification:

| Sum of points | Category |
|---------------|----------|
| >= 10 | **Pathogenic** |
| 6-9 | **Likely Pathogenic** |
| 0-5 | **VUS** |
| -1 to -6 | **Likely Benign** |
| <= -7 | **Benign** |

InterVar / GeneBe / VarSome / Franklin all implement Tavtigian point summation under the hood. Combinations never appearing in the 2015 combining rules (e.g., PVS1_VeryStrong + PM2_Supporting -> LP) emerge naturally from point arithmetic.

## PVS1 Decision Tree (Abou Tayoun 2018 *Hum Mutat* 39:1517)

PVS1 is the most consequential code: pathogenic Very Strong (8 points) for predicted loss-of-function in a gene where LoF is established disease mechanism. The 2018 decision tree refined PVS1 from a binary into a graded code based on:

1. **Variant type**: nonsense / frameshift / canonical +-1,2 splice / initiation codon / single-exon deletion / multi-exon deletion.
2. **NMD prediction**: variant in 5'-most exon OR >50bp upstream of last exon-exon junction -> NMD-triggered. Else truncated protein.
3. **Critical region**: removal of >10% of coding sequence OR removal of a critical functional domain.
4. **Alternative isoform**: does the variant affect a transcript expressed in disease-relevant tissue?

Output strengths:

| Output | Original Strength |
|--------|-------------------|
| PVS1_VeryStrong | Strongest (Very Strong) |
| PVS1_Strong | Strong |
| PVS1_Moderate | Moderate |
| PVS1_Supporting | Supporting |

**Subsumption rule** (Abou Tayoun 2018): PVS1 + PP3 -> only PVS1 counts (PP3 is subsumed). Same for PVS1 + PM4.

**>15 VCEP-specific PVS1 trees exist** as of 2024 (CDH1, ENIGMA BRCA1/2, FH LDLR/APOB/PCSK9, InSiGHT MMR, RASopathies, hearing loss, hypertrophic cardiomyopathy, Rett/Angelman, etc.). The automated implementation is **AutoPVS1** (Xiang 2020).

## Pejaver 2022 PP3/BP4 Calibrated Thresholds (the load-bearing 2024+ calibration)

Pejaver 2022 *AJHG* 109:2163 Bayesian-calibrated 13 missense predictors to PP3/BP4 strength levels using ClinVar P/B variants with leave-one-gene-out cross-validation.

| Predictor | BP4_Strong | BP4_Moderate | BP4_Supporting | PP3_Supporting | PP3_Moderate | PP3_Strong | Fails when |
|-----------|-----------|--------------|----------------|----------------|--------------|------------|-----------|
| **REVEL** | <= 0.016 | <= 0.183 | <= 0.290 | >= 0.644 | >= 0.773 | >= 0.932 | Stacked with BayesDel/VEST4 (training overlap; double-counting) |
| **BayesDel (no AF)** | n/a | <= -0.36 | <= -0.18 | >= 0.13 | >= 0.27 | >= 0.50 | No BP4_Strong reached; combine no-AF version with PM2_Supporting |
| **VEST4** | n/a | <= 0.302 | <= 0.449 | >= 0.764 | >= 0.861 | >= 0.965 | No BP4_Strong reached; indels (missense-trained); regulatory variants |
| **MutPred2** | (Pejaver 2022) | -- | -- | -- | -- | -- | Genes with sparse MAVE training data |
| **AlphaMissense** | no -4 interval | <= 0.099 (with a separate -3 interval <= 0.070) | 0.100-0.169 | 0.792-0.905 | 0.906-0.971 | >= 0.990 (with a separate +3 interval 0.972-0.989) | Developer threshold 0.564 misapplied as PP3 |

**The numbers to memorize: REVEL >= 0.932 = PP3_Strong; REVEL <= 0.016 = BP4_Strong (<= 0.003 BP4_VeryStrong); the 0.290-0.644 band is indeterminate (no criterion).**

**AlphaMissense calibration** (Bergquist 2025 *Genet Med* 27:101402): `<=0.070` is -3 points, `0.071-0.099` BP4 Moderate, `0.100-0.169` BP4 Supporting, `0.170-0.791` indeterminate, `0.792-0.905` PP3 Supporting, `0.906-0.971` PP3 Moderate, `0.972-0.989` +3 points, and `>=0.990` PP3 Strong. The +/-3 bands belong to the published point-based table and must not be rounded to a legacy strength. The developer threshold 0.564 is inside the calibrated indeterminate interval.

**Do not stack predictors.** REVEL, BayesDel, VEST4 share ClinVar/HGMD training data; using REVEL >= 0.773 AND BayesDel >= 0.27 to claim "two independent moderate hits" is double-counting. Pejaver 2022 explicitly recommends using ONE predictor per variant.

## PM2_Supporting (ClinGen SVI 2020)

The original PM2 ("absent from controls") was over-weighted. SVI 2020 downgraded to **PM2_Supporting** (1 point, not 2). Mechanism: most rare variants are benign. Empirical recalibration showed ~6 variants per gene downgrade from LP to VUS when PM2 -> Supporting. Many 2017-2019 LP curations require re-classification post-SVI 2020 update.

## PS3/BS3 Functional Evidence (Brnich 2020 *Genome Med* 12:3)

OddsPath framework; the four-step SOP:

1. Define disease mechanism for the gene.
2. Evaluate assay class (e.g., MAVE, biochemical, animal model).
3. Evaluate specific assay instance (controls, replicate consistency).
4. Apply per-variant.

OddsPath calibration mapping to ACMG strengths:

| OddsPath | Pathogenic strength | Benign strength |
|----------|--------------------:|----------------:|
| > 350 | Very Strong | n/a |
| > 18.7 to 350 | Strong | < 0.053 = Strong |
| > 4.3 to 18.7 | Moderate | 0.053 to < 0.23 = Moderate |
| > 2.1 to 4.3 | Supporting | 0.23 to < 0.48 = Supporting |
| 0.48 to 2.1 | Indeterminate | Indeterminate |

**MAVEdb deep-mutational scans** with >=11 controls (>=5 P/LP + >=5 B/LB) can yield up to PS3_Strong/BS3_Strong via OddsPath calibration. This is the entry point for MAVE/saturation-mutagenesis evidence into ACMG.

**Default-Strong PS3 application is over-strengthening** without OddsPath calibration. Under Brnich Table 3, OddsPath must exceed 18.7 for Strong and exceed 350 for Very Strong; an uncalibrated assay needs assay-validity review rather than a guessed strength.

## ClinGen SVI Splicing Subgroup 2023 (Walker *AJHG* 110:1046)

**SpliceAI is the recommended primary splicing tool.** Calibrated thresholds:

| SpliceAI DS_max | Strength (Walker 2023: computational splice codes applied at Supporting weight) |
|-----------------|----------|
| >= 0.2 | PP3_Supporting (minimum threshold for ANY splicing PP3) |
| 0.1 - 0.2 | Indeterminate (no criterion) |
| <= 0.1 | BP4_Supporting |

SpliceAI prediction alone does NOT reach PP3_Strong; strength escalation requires experimental/RNA splicing evidence (PS3) or the repurposed PVS1 route.

**SpliceVault / 300K-RNA** (Dawes 2023 *Nat Genet* 55:324): does NOT predict whether a variant is splice-altering; predicts WHAT the aberrant transcript will be (which exon skips, which cryptic site activates). 96% sensitivity for exon-skipping; 86% for cryptic site activation in 140 clinical RNA-tested cases. Critical for PVS1 application to splice variants because PVS1 depends on whether the aberrant transcript triggers NMD.

**Pangolin** (Zeng 2022 *Genome Biol* 23:103): SpliceAI improvement for cryptic donor sites; not yet ClinGen-endorsed but increasingly used as tiebreaker.

## BS1 / BA1 (Whiffin Max-Credible-AF)

BA1 default: AF > 5% in non-bottleneck group per ClinGen SVI; VCEP-specific overrides (Hearing Loss VCEP uses 0.5% AR).

BS1 gene-specific: `(prevalence x heterogeneity x allelic-contribution) / (penetrance x 2)` from Whiffin 2017 *Genet Med* 19:1151. Compare against gnomAD `grpmax_faf95`.

See `clinical-databases/gnomad-frequencies` for FAF95 details.

## ClinGen VCEP CSpec Hierarchy

| Layer | Authority | Application |
|-------|-----------|-------------|
| Generic ACMG/AMP 2015 | Richards 2015 | Default fallback |
| ClinGen SVI specifications | SVI Working Group | Overrides generic for all genes (PM2 -> Supporting; AutoPVS1 trees; etc.) |
| VCEP-specific CSpec | Gene/disease-specific expert panel | Overrides SVI for that gene-disease |

**ClinGen VCEP CSpec current REST workflow:** inspect `https://cspec.genome.network/cspec/srvc`, resolve the gene identity, list `/cspec/Gene/id/{GENE}/SequenceVariantInterpretation/version`, and retain the exact version/access date. Do not use the retired `/cspec/ui/svi/all` GET route. Examples of gene-specific overrides include:
- Hearing Loss VCEP: PM2 -> supporting default; PS3 thresholds upgraded for OTOF; BA1 lowered to 0.5% AR.
- ENIGMA BRCA1/2 VCEP: gene-specific PVS1 trees with NMD escape rules; PS4 case-control thresholds.
- Inherited Cardiac Conditions VCEP: gene-specific PS4 (5+ unrelated probands for PS4_Supporting).

**Apply VCEP CSpec when one exists.** Generic ACMG with no VCEP awareness is unreliable for many genes.

## Cancer Somatic Framework (Li 2017 *J Mol Diagn* 19:4)

AMP/ASCO/CAP somatic variant interpretation; four tiers:

| Tier | Definition | Action |
|------|-----------|--------|
| **Tier I-A** | FDA-approved drug for same tumor type with this biomarker | On-label therapy |
| **Tier I-B** | Professional guidelines (NCCN, ESMO) | Standard-of-care |
| **Tier II-C** | FDA drug in different tumor type (off-label) | Basket trials |
| **Tier II-D** | Preclinical / investigational | Research |
| **Tier III** | VUS-somatic | Watch list |
| **Tier IV** | Benign-somatic | Filter out |

**Knowledgebases:** OncoKB (MSKCC; Chakravarty 2017), CIViC (Griffith 2017 *Nat Genet* 49:170), CGI (Tamborero 2018), JAX-CKB, COSMIC. Knowledgebase levels do not map one-to-one to AMP tiers. Record the evidence record and access date, distinguish regulatory/guideline, clinical, and preclinical evidence, and separately record clinical significance. Tier III means uncertain significance; Tier IV means likely benign/benign. Do not infer a tier from a bare numeric OncoKB level.

The **Variant Interpretation for Cancer Consortium (VICC) Meta-Knowledgebase** standards (2024-2025) harmonize across knowledgebases. ClinGen Somatic VCEPs are emerging (started 2022).

## Decision Tree by Variant Type

| Variant type | Recommended workflow |
|--------------|----------------------|
| Predicted LoF in known LoF-mechanism gene | AutoPVS1 decision tree -> PVS1_VeryStrong/Strong/Moderate/Supporting; check VCEP-specific PVS1 |
| Missense in known missense-pathogenic gene | Apply Pejaver 2022 PP3/BP4 calibrated thresholds; ONE predictor only |
| Splice variant | SpliceAI DS_max + SpliceVault for aberrant-transcript prediction; PP3_Supporting if DS_max >=0.2 (strength escalation needs RNA/experimental evidence) |
| Synonymous | SpliceAI for cryptic splice effect; synVep / PrimateAI synonymous extension |
| Variant in ACMG SF v3.2 gene | Apply full classification; flag P/LP for opt-in disclosure |
| Cancer somatic variant | AMP/ASCO/CAP 2017 Tier I-IV; cross-check OncoKB / CIViC |
| Variant in Limited gene-disease validity | ClinGen Strong/Definitive required for clinical action |
| Functional evidence available | Brnich 2020 PS3/BS3 OddsPath framework |
| Family segregation | PP1 / BS4 LOD score per Biesecker 2024 |
| In-trans observations (AR) | PM3 with ClinGen tabular scoring system |
| HGVS-c on alternative transcript | Re-evaluate on MANE Select |

## Standard Workflow: ACMG Classification

**Goal:** Apply ACMG/AMP framework to a candidate variant with proper SVI specifications and VCEP overrides.

**Approach:** Pull aggregated evidence; apply Pejaver-calibrated in-silico thresholds; check VCEP-specific CSpec; sum Tavtigian points.

Run the reusable implementation in [`../scripts/acmg_classify.py`](../scripts/acmg_classify.py). It retains the source example's Tavtigian scoring, calibrated predictor helpers, PVS1 sketch, OddsPath helper, population-frequency helper, GeneBe request, and somatic-tier helper.

## Per-Operation Failure Modes

**1. Stacking REVEL + BayesDel + VEST4 as independent evidence**
- Trigger: Apply PP3 from multiple predictors.
- Mechanism: Predictors share training data; double-counting.
- Symptom: Inflated PP3 strength; over-classified LP.
- Fix: Use ONE predictor per variant (Pejaver 2022).

**2. AlphaMissense evidence with developer 0.564 threshold**
- Trigger: Apply AlphaMissense >0.564 -> PP3_Strong.
- Mechanism: 0.564 is the developer-recommended likely-pathogenic threshold, NOT a calibrated PP3 cutoff; use the ClinGen SVI calibration (Bergquist 2025) score thresholds instead.
- Symptom: Over-application of PP3.
- Fix: Use Bergquist 2025 calibrated intervals; 0.564 is indeterminate, Supporting begins at 0.792, Moderate at 0.906, and Strong at 0.990.

**3. PVS1 applied to nonsense variant in GoF gene**
- Trigger: Nonsense variant in SCN5A reported as PVS1 for LQT3.
- Mechanism: SCN5A has both LoF (Brugada) and GoF (LQT3) mechanisms. PVS1 should NOT apply if LoF is not the established mechanism.
- Symptom: Wrong classification; clinical action mis-directed.
- Fix: Check ClinGen gene-disease mechanism; apply PVS1 only when LoF is established.

**4. Generic ACMG instead of VCEP CSpec**
- Trigger: Apply default ACMG to a variant in a gene with VCEP-specific CSpec.
- Mechanism: VCEP CSpec overrides for the gene; e.g., Hearing Loss VCEP PM2 default = supporting, BA1 = 0.5%.
- Symptom: Wrong strength applied; misclassification.
- Fix: Use the current CSpec REST service, resolve the gene's version list, and record the selected gene-specific specification plus access date.

**5. PM2 at Moderate (pre-2020 SVI)**
- Trigger: Apply PM2 = Moderate (2 points) per 2015 rules.
- Mechanism: SVI 2020 downgraded to PM2_Supporting (1 point).
- Symptom: ~6 variants per gene over-strengthened LP.
- Fix: Use PM2_Supporting per current SVI.

**6. PS3 default Strong without OddsPath**
- Trigger: Apply PS3 = Strong without OddsPath calibration.
- Mechanism: Default PS3 = Strong over-strengthens; Brnich 2020 SOP requires OddsPath > 4.3 for Moderate and > 18.7 for Strong.
- Symptom: PS3-driven over-classification.
- Fix: Apply Brnich 2020 four-step OddsPath; >4.3 is Moderate, >18.7 is Strong, and >350 is Very Strong. Without calibration, stop for assay-validity review.

**7. Synonymous treated as no impact**
- Trigger: Filter out synonymous variants from classification pipeline.
- Mechanism: Synonymous can disrupt splicing; SpliceAI captures this.
- Symptom: Pathogenic splice-disrupting synonymous missed.
- Fix: Always run SpliceAI on synonymous variants in disease genes; PP3_Supporting if DS_max >= 0.2.

**8. ClinVar P + ClinGen Limited validity**
- Trigger: Report variant P in gene with Limited gene-disease validity.
- Mechanism: ClinVar P is variant-level; gene-disease validity is upstream.
- Symptom: Mis-attribution to a non-disease gene.
- Fix: Apply ClinGen gene-disease validity gate (Moderate+ for clinical action); for Limited genes, require VCEP curation.

**9. Variant on wrong transcript**
- Trigger: HGVS-c on alt transcript; functional impact different on MANE Select.
- Mechanism: Tissue-specific isoform considerations; MANE Select 2024+ is clinical standard.
- Symptom: Wrong consequence prediction.
- Fix: Re-evaluate on MANE Select transcript; cross-check with VEP `--mane_select`.

## Reconciliation: When Tools Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| GeneBe LP vs VarSome P | Different VCEP-specific application | Check VCEP CSpec; apply gene-specific rules |
| ClinVar P vs my classification VUS | Submission stale OR my evidence incomplete | Re-curate with current evidence; check ClinVar star + freshness |
| REVEL PP3_Strong vs SpliceAI BP4 | Variant has missense impact but no splice impact | Apply ONE predictor; if splice-altering, PVS1 trumps |
| PVS1 applies but ClinGen Limited validity | Variant-level vs gene-disease tension | Treat as candidate; require VCEP or strong functional evidence |
| ClinGen VCI vs automated tool | VCI is gold standard for expert curation | Trust VCI; automated tools approximate |
| AlphaMissense 0.564 dev call vs calibrated strength | Developer threshold lies in the calibrated indeterminate interval | Apply Bergquist 2025 exact intervals |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| Tavtigian P | >= 10 points | Tavtigian 2020 |
| Tavtigian LP | 6-9 points | Tavtigian 2020 |
| Tavtigian VUS | 0-5 points | Tavtigian 2020 |
| Tavtigian LB | -1 to -6 | Tavtigian 2020 |
| Tavtigian B | <= -7 | Tavtigian 2020 |
| REVEL PP3_Strong | >= 0.932 | Pejaver 2022 |
| REVEL BP4_Strong | <= 0.016 | Pejaver 2022 |
| SpliceAI PP3 (Supporting) | >= 0.2 | Walker 2023 |
| SpliceAI BP4 (Supporting) | <= 0.1 | Walker 2023 |
| BA1 default | grpmax_faf95 > 5% | ClinGen SVI |
| BS1 | grpmax_faf95 > gene-specific max-credible-AF | Whiffin 2017 |
| PM2 -> PM2_Supporting | Always (post-SVI 2020) | SVI 2020 |
| PS3 OddsPath Supporting / Moderate / Strong / Very Strong | >2.1 / >4.3 / >18.7 / >350 | Brnich 2020 |
| PVS1 LoF mechanism check | Required (do not apply if GoF) | Abou Tayoun 2018 |
| ACMG SF v3.2 | 81 genes | Miller 2023 |
| Cancer Tier I-A | FDA drug + same tumor + this biomarker | Li 2017 |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| Over-application of PP3 | Multiple predictors stacked | ONE predictor only |
| AlphaMissense PP3_Strong from dev threshold | 0.564 not calibrated | Use Pejaver-style REVEL |
| LP variant in gene with Limited validity | No gene-disease gate | Apply ClinGen gene-disease validity |
| PVS1 in GoF gene | Wrong mechanism | Check ClinGen gene-disease mechanism |
| Non-VCEP rule for VCEP-covered gene | Generic ACMG | Apply VCEP CSpec |
| PM2 = Moderate | Pre-SVI 2020 | Use PM2_Supporting |
| PS3 = Strong at OddsPath 8 | Shifted thresholds | OddsPath 8 is Moderate; Strong begins above 18.7 |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why Tavtigian point system?" | Every modern automated classifier implements it (InterVar, GeneBe, VarSome, Franklin). The 2015 combining rules are subsumed; many P/LP combinations only emerge from points. |
| "Why ONE predictor and not REVEL + BayesDel?" | Pejaver 2022 explicit recommendation; predictors share training data. |
| "AlphaMissense PP3_Strong?" | Bergquist 2025: Strong begins at 0.990; 0.972-0.989 is +3 points and 0.564 is indeterminate. |
| "PVS1 for nonsense in SCN5A LQT3" | LQT3 is GoF; LoF mechanism not established; PVS1 does not apply. |
| "Generic ACMG vs VCEP" | VCEP CSpec overrides generic; we check `cspec.genome.network` for active VCEP. |
| "Splice variant PP3 from SpliceAI" | Walker 2023 SVI Splicing Subgroup: DS_max >= 0.2 applies PP3 at Supporting weight; prediction alone does not reach PP3_Strong (needs RNA/experimental evidence). |
| "PM2 Moderate or Supporting?" | SVI 2020 downgraded to Supporting; we use Supporting for all classification post-2020. |

## References

- Richards S et al. 2015. Standards and guidelines for the interpretation of sequence variants. *Genet Med* 17:405. (Original ACMG/AMP)
- Tavtigian SV et al. 2018. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med* 20:1054.
- Tavtigian SV et al. 2020. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Hum Mutat* 41:1734.
- Abou Tayoun AN et al. 2018. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. *Hum Mutat* 39:1517.
- Pejaver V et al. 2022. Calibration of computational tools for missense variant pathogenicity classification. *Am J Hum Genet* 109:2163.
- Bergquist T et al. 2025. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 27:101402.
- Brnich SE et al. 2020. Recommendations for application of the functional evidence PS3/BS3 criterion using the ACMG/AMP sequence variant interpretation framework. *Genome Med* 12:3.
- Walker LC et al. 2023. ClinGen SVI Splicing Subgroup recommendations. *Am J Hum Genet* 110:1046.
- Biesecker LG et al. 2024. ClinGen guidance for use of the PP1/BS4 co-segregation and PP4 phenotype specificity criteria for sequence variant pathogenicity classification. *Am J Hum Genet* 111:24. (PP1/BS4 co-segregation)
- Cheng J et al. 2023. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 381:eadg7492.
- Jaganathan K et al. 2019. Predicting splicing from primary sequence with deep learning. *Cell* 176:535. (SpliceAI)
- Zeng T et al. 2022. Predicting RNA splicing from DNA sequence using Pangolin. *Genome Biol* 23:103.
- Dawes R et al. 2023. SpliceVault predicts the precise nature of variant-associated mis-splicing. *Nat Genet* 55:324.
- Whiffin N et al. 2017. Using high-resolution variant frequencies to empower clinical genome interpretation. *Genet Med* 19:1151.
- Li MM et al. 2017. Standards and guidelines for the interpretation and reporting of sequence variants in cancer. *J Mol Diagn* 19:4. (AMP/ASCO/CAP)
- Miller DT et al. 2023. ACMG SF v3.2 list. *Genet Med* 25:100866.
- Stawiński P, Płoski R. 2024. Genebe.net: implementation and validation of an automatic ACMG variant pathogenicity criteria assignment. *Clin Genet* 106:119.
- Kopanos C et al. 2019. VarSome: the human genomic variant search engine. *Bioinformatics* 35:1978.
- Li Q, Wang K. 2017. InterVar: clinical interpretation of genetic variants. *Am J Hum Genet* 100:267.
- Xiang J et al. 2020. AutoPVS1 -- automated PVS1 decision-tree implementation (verify exact venue/year against the published code/release).
- ClinGen CSpec Registry REST service: `https://cspec.genome.network/cspec/srvc`
- ClinGen VCI (Variant Curation Interface): `https://curation.clinicalgenome.org/`

