---
name: bio-causal-genomics-effector-gene-prioritization
description: Maps GWAS-implicated loci to candidate effector (causal) genes by integrating variant-to-gene (V2G) features via Open Targets L2G (Mountjoy 2021), MAGMA gene-based association (de Leeuw 2015), FUMA SNP2GENE, cS2G combined SNP-to-gene scores (Gazal 2022), Polygenic Priority Scores (PoPS, Weeks 2023), FLAMES, INQUISIT, DEPICT, and enhancer-gene predictors (ABC, ENCODE-rE2G). Use when narrowing a GWAS lead locus to a candidate causal gene, picking between proximity, eQTL-based, and similarity-based prioritizers, integrating multi-evidence streams (fine-mapping, colocalization, ABC enhancer-gene, distance, chromatin), reconciling discordant L2G vs PoPS calls, prioritizing tissue-specific eQTL evidence, or triangulating across at least three independent lines of evidence for a publication-grade effector-gene nomination.
tool_type: mixed
primary_tool: MAGMA
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: MAGMA 1.10+ (cncr.nl/research/magma), FUMA web platform v1.6+ (fuma.ctglab.nl), Open Targets Genetics API (REST + GraphQL, June 2024 release), PoPS (head of `FinucaneLab/pops`, 2024), cS2G pre-computed scores (Zenodo record 7754032, Gazal 2022), ABC-Enhancer-Gene-Prediction 0.2.2+, ENCODE-rE2G v1.0+ (Gschwind 2023 preprint), DEPICT v1 rel194, INQUISIT (Fachal 2020 supplementary), Python 3.9-3.11, R 4.3+, PLINK 1.9 + PLINK 2.0.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `magma --help` to confirm gene-window, gene-annot, and gene-set flag names
- Python: `pip show gentropy`; introspect endpoints at `api.platform.opentargets.org/api/v4/graphql`
- R: `packageVersion('coloc')` etc. for upstream evidence integration

If a script throws an error about an argument that has moved (e.g. an Open Targets endpoint renamed during a release) or a model file schema change, introspect the installed tool and adapt rather than retrying. Open Targets Genetics deprecated the standalone Genetics Portal in 2024 in favour of the integrated platform; verify endpoint URLs at the time of use.

# Effector Gene Prioritization

**"Which gene at this GWAS locus is actually the causal mediator?"** -> Integrate fine-mapping, colocalization, chromatin-based enhancer-gene predictions, distance, and gene-similarity priors into a per-locus per-gene confidence score, then require concordance across multiple orthogonal evidence streams before nominating a causal effector. Effector gene prioritization is the bridge between statistical fine-mapping (variant level) and biological hypothesis (gene level); it is the most failure-prone step in GWAS-to-target pipelines because the nearest-gene assumption is wrong roughly 30-50% of the time at well-studied loci.

- CLI (gene-level association): `magma --bfile ref --gene-loc geneloc.txt --pval gwas.tsv ncol=N --out out` -> `magma --gene-results out.genes.raw --set-annot annot.txt --out out`
- Web (integrative): FUMA SNP2GENE at fuma.ctglab.nl (positional + eQTL + Hi-C + chromatin in one workflow)
- API (pre-computed L2G): Open Targets Genetics GraphQL `studyLocus2GeneTable` query (note: Open Targets Genetics was consolidated into the Open Targets Platform in 2024; verify the live endpoint at `api.platform.opentargets.org/api/v4/graphql`)
- Python (similarity prior): `python pops.py --gene_annot_path gene_annot.txt --feature_mat_prefix features --control_features_path control.features --magma_prefix magma_out --out_prefix out`
- Lookup (combined SNP-to-gene): cS2G pre-computed gene scores at zenodo.org/records/7754032
- CLI (enhancer-gene): ABC pipeline or ENCODE-rE2G (cross-reference atac-seq/enhancer-gene-linking)

V2G is not one method but a portfolio. Open Targets L2G aggregates per-locus per-gene features (distance + coloc + chromatin + V2G) trained on curated gold-standard genes; PoPS adds an orthogonal genome-wide polygenic prior from gene-pathway co-membership; MAGMA provides the lightweight gene-level p-value baseline. Strong effector calls emerge from concordance across these orthogonal signal types, not from any single tool.

## Scope

Every method here nominates a candidate causal gene at the population level, from GWAS
summary statistics -- it is hypothesis-generating research output, not an individual clinical
determination. A gene scoring as a locus's effector (even a well-known drug-target gene like
PCSK9, LDLR, or APOE) says nothing about which treatment is right for a specific patient with a
specific variant, comorbidities, and clinical history. If a request asks this Skill's output to
drive a per-patient prescribing, diagnostic, or treatment-selection decision, decline that framing
explicitly and redirect to the patient's treating clinician or a clinical genetics / pharmacogenomics
pathway; offer the population-level evidence (concordance tier, evidence streams) as context for that
clinician, not as a substitute for their judgment.

## Prerequisites

- GWAS summary statistics (SNP, CHR, BP, A1, A2, BETA or Z, SE, P, N)
- Fine-mapping output (SuSiE / FINEMAP credible sets) from causal-genomics/fine-mapping
- Colocalization output (coloc.abf or coloc.susie PP.H4) from causal-genomics/colocalization-analysis
- Ancestry-matched LD reference panel (e.g. 1000 Genomes EUR PLINK bfile)
- Pre-computed integrative scores when available: Open Targets L2G, cS2G, PoPS feature matrix
- For matched-tissue enhancer-gene linking: ATAC-seq + H3K27ac ChIP + Hi-C/Micro-C in the candidate causal tissue (cross-reference atac-seq/enhancer-gene-linking)

### eQTL / pQTL Panel Selection

Panel choice dominates coloc and TWAS sensitivity (see "eQTL tissue mis-specification" in `references/failure-modes.md`). Match to the causal tissue identified via LDSC-SEG or stratified LDSC before locking in.

| Panel | N donors | Tissues / cells | Use case |
|-------|----------|-----------------|----------|
| GTEx v8 | 838 | 49 tissues | PredictDB standard for 2026 TWAS pipelines |
| GTEx v10 (2024) | ~1,000 | 49 tissues | Limited PredictDB coverage; transitional |
| eQTLGen | 31,684 | Whole blood | Maximum statistical power for blood-relevant traits |
| eQTLGen Phase 2 (sc) | ~30,000 | sc immune | Cell-type-resolved blood eQTLs |
| OneK1K | 982 | PBMC, 14 cell types | sc-eQTL discovery for immune traits |
| INTERVAL pQTL | 3,301 | Plasma proteome | Drug-target cross-reference via pQTL-MR |

## Algorithmic Taxonomy

Ten tools compete for this job: Open Targets L2G, V2G, MAGMA, FUMA SNP2GENE, cS2G, PoPS, FLAMES,
INQUISIT, DEPICT, and ABC/ENCODE-rE2G (plus sc-eQTL/cell-type TWAS). Each has a distinct model,
input/output shape, and failure regime. See `references/algorithmic-taxonomy.md` for the full
per-tool comparison table (model, inputs, output, strength, fails-when) before picking or
contrasting specific tools. The load-bearing summary: L2G + PoPS is the current de facto
two-method baseline (orthogonal by construction -- see Reconciliation below); cS2G is the
heritability-calibrated per-SNP lookup; FUMA is the no-install community standard. Methodology
evolves; verify against the current Open Targets release (platform-docs.opentargets.org), the
latest PoPS feature matrix at FinucaneLab/pops, and ABC / ENCODE-rE2G releases before locking
on a single prioritiser.

## Decision Tree by Scenario

| Scenario | Recommended workflow | Why |
|----------|---------------------|-----|
| Open Targets Platform covers the trait | Query L2G via GraphQL + cross-check V2G; sanity-check with PoPS | Pre-computed, gold-standard-validated; minimal compute |
| Custom trait, EUR GWAS sumstats only | MAGMA + manual fine-mapping (SuSiE) + coloc per QTL panel | Build evidence streams from primitives; combine in own integrative scorer |
| Tissue known (e.g. liver for lipid traits) | Tissue-specific eQTL coloc + ABC / ENCODE-rE2G + S-PrediXcan + L2G | Tissue-targeted evidence reduces false positives from wrong-tissue eQTLs |
| Tissue unknown a priori | LDSC-SEG (Finucane 2018 Nat Genet 50:621) to prioritise tissue + S-MultiXcan + PoPS | Identify causal tissue before locking on a single eQTL panel |
| Distal / long-range regulation suspected | ABC / ENCODE-rE2G / HiChIP / Cicero overlay; deprioritise distance-only methods | Nearest-gene fails ~ 30-50% of the time at well-studied loci |
| Polygenic background trait (e.g. height, BMI) | L2G + PoPS concordance | PoPS captures pathway prior absent from L2G features; concordance flags strong candidates |
| Publication-grade triangulation | All evidence streams; require concordance across >= 3 orthogonal lines | High-confidence claim defensible to reviewers |
| Multi-ancestry GWAS | MAGMA per ancestry + ancestry-specific eQTL coloc + MA-FOCUS for TWAS | Single-ancestry weights miscalibrated for other ancestries |
| Locus with no coding variants, no significant eQTL | ABC / ENCODE-rE2G in candidate tissue + chromatin annotation; tag as "regulatory of unknown gene" | Distance + chromatin may be the only signal; acknowledge low confidence |
| HLA region (chr6:28477797-33448354 hg19; chr6:28510120-33480577 hg38; extended chr6:25-35 Mb both builds) | Exclude or use HLA-imputation; do not run standard V2G; verify build before excluding | Long-range LD breaks every gene-by-gene method |

## Per-Method Failure Modes

Six recurring pitfalls, each with trigger / mechanism / symptom / fix: the nearest-gene
assumption (wrong 30-50% of the time), eQTL tissue mis-specification, MAGMA gene-window choice,
coloc failing at multi-causal-variant loci, PoPS-vs-L2G discordance, and pleiotropic loci with
multiple true effector genes. See `references/failure-modes.md` for the full write-up of each --
consult it when diagnosing a specific discordance rather than before every run.

## Per-Credible-Set Gene-Assignment Hierarchy

When fine-mapping returns credible sets (cross-reference causal-genomics/fine-mapping), each credible variant should be assigned to a gene using a fixed lexicographic ladder rather than a single feature. The ladder collapses ambiguity by preferring direct mechanistic evidence first and falling back to weaker signals only when stronger ones are absent:

1. **Coding consequence at credible variant** (missense, splice-donor / splice-acceptor, stop-gained, start-lost via VEP / Ensembl consequence) -> assign variant to that gene.
2. **eQTL / pQTL colocalization PP.H4 >= 0.7** with the gene's expression QTL (matched tissue) -> assign to that gene (cross-reference causal-genomics/colocalization-analysis).
3. **ABC or ENCODE-rE2G enhancer-gene linkage** when matched ATAC + H3K27ac (+ optional Hi-C) is available in the candidate tissue -> assign to the linked gene (cross-reference atac-seq/enhancer-gene-linking).
4. **Nearest TSS** -> assign as last-resort fallback; flag as low-confidence (nearest-gene is correct only 50-70% of the time at well-fine-mapped loci).

If a single credible variant ties across two or more genes at the same rung (e.g. coding consequence in gene A AND a competing eQTL coloc to gene B), report multi-gene candidacy explicitly rather than forcing a single assignment; the locus may be multi-effector or the credible set may straddle a regulatory boundary.

## Multi-Evidence Integration Framework

A strong candidate causal gene at a GWAS locus requires concordance across multiple orthogonal evidence streams. The six canonical streams:

| Evidence stream | What it tests | Pass threshold | Source |
|----------------|---------------|----------------|--------|
| Fine-mapping | SuSiE PIP for variant; variant assigned to gene by ABC or eQTL coloc | PIP > 0.5; credible set purity > 0.5 | susieR (cross-reference causal-genomics/fine-mapping) |
| Colocalization | Shared causal variant with gene's eQTL / pQTL / sQTL | coloc.abf or coloc.susie PP.H4 >= 0.7 | coloc (cross-reference causal-genomics/colocalization-analysis) |
| Distance | Variant within annotated regulatory unit (gene body or enhancer-gene unit) | Distance to TSS <= 100 kb OR within ABC enhancer-gene unit | Convention; Mountjoy 2021 |
| Polygenic prior (similarity) | Gene shares pathway / co-expression / PPI with other GWAS hits for the trait | PoPS score in top decile per locus | Weeks 2023 |
| L2G (integrative classifier) | Per-locus per-gene gradient-boosting on a panel of features | L2G score >= 0.5 (Open Targets default high-confidence) | Mountjoy 2021 |
| Chromatin / enhancer-gene | ABC or ENCODE-rE2G connects fine-mapped variant to gene's promoter | ABC >= 0.02 OR ENCODE-rE2G >= 0.5 | Fulco 2019; Gschwind 2023 |
| Deep-learning variant effect (optional 7th) | chromBPNet / EnFormer in silico variant effect on accessibility/expression at credible variant | \|log2FC\| > 1 (chromBPNet strong-effect); agreement across two models | chromBPNet: Pampari 2025 bioRxiv 2024.12.25.630221 (preprint); EnFormer: Avsec 2021 Nat Methods 18:1196; cross-reference atac-seq/deep-learning-atac |
| Cross-trait coincidence (optional 8th) | Same gene flagged at related-trait loci (e.g. lipid GWAS at CHD lead variant) | Same gene top-ranked at >= 2 related traits | Convention |

**Operational rule:** Report a gene as a "high-confidence causal effector" only when >= 3 of the 6 core evidence streams are concordant (deep-learning variant effect and cross-trait coincidence are optional supplementary streams). >= 4 concordant is "strong-confidence"; >= 5 is "near-certain pending experimental validation". Single-stream evidence is associational only; two-stream concordance is suggestive. This is the standard used by Open Targets Genetics (Mountjoy 2021), GTEx-derived target nomination pipelines (Open Targets Platform), and pharma drug-discovery workflows.

## Quantitative Thresholds

Canonical pass/fail thresholds for every evidence stream (L2G, MAGMA gene-wide p, coloc PP.H4,
PoPS decile, ABC/ENCODE-rE2G, cS2G, distance-to-TSS, MAGMA window, fine-mapping PIP,
multi-evidence concordance, single-cell eQTL panel size) live in
`references/quantitative-thresholds.md`. Consult it when scoring or reporting a specific
evidence stream -- the two thresholds used constantly are repeated inline where needed
(coloc PP.H4 >= 0.7 in the Multi-Evidence table below; >= 3 of 6 concordance in the
Operational rule).

## MAGMA Gene-Based and Gene-Set Pipeline

**Goal:** Compute gene-level p-values from GWAS summary statistics and test gene sets (e.g. MSigDB pathways) for enrichment.

**Approach:** Pre-format the SNP-to-gene annotation (per-gene SNP membership using a configurable window); run gene-based analysis with `--gene-results`; downstream, test gene sets via `--set-annot`. MAGMA's lambda-correction handles LD via the reference panel.

```bash
# Step 1: SNP-to-gene annotation using a 35kb upstream + 10kb downstream window (FUMA default)
magma --annotate window=35,10 \
    --snp-loc gwas.snploc \
    --gene-loc NCBI37.3.gene.loc \
    --out annot_35_10

# Step 2: Gene-based association (raw GWAS sumstats; multi-model approach)
magma --bfile g1000_eur \
    --pval gwas.pval ncol=N \
    --gene-annot annot_35_10.genes.annot \
    --out gene_step

# Step 3: Gene-set enrichment via competitive testing (recommended over self-contained)
magma --gene-results gene_step.genes.raw \
    --set-annot msigdb_v7.5_C2.gmt \
    --out gene_set_step

# Inspect: gene_step.genes.out (per-gene Z, p); gene_set_step.gsa.out (per-set p)
```

The `--gene-annot` window choice is the dominant methodological lever; 35kb upstream + 10kb downstream is the FUMA recommendation but is not universally accepted. Sensitivity over 0+0, 35+10, and 50+50 is good practice for high-stakes reports.

**Step 3 has an undocumented minimum gene count.** MAGMA's `--set-annot` competitive regression conditions on 6 internal covariates (gene size, log(gene size), gene density, log(gene density), inverse MAC, log(inverse MAC)). With too few genes in `--gene-results` relative to those 6 covariates, the regression is structurally unidentifiable: MAGMA aborts with `ERROR: insufficient degrees of freedom to run analyses` (verified on a real 5-gene locus run), or, with even fewer genes, fails earlier with `no input variables to analyse` from a gene-set-variance check (verified on a real 3-gene run). Rule of thumb: gene-set enrichment needs several hundred+ genes -- **a single-locus MAGMA run (a handful of genes) should skip Step 3 entirely**; gene-set enrichment is a genome-wide- or many-loci-scale analysis, not a per-locus one. `examples/magma_genebased.sh` checks the gene count and skips Step 3 automatically below this threshold. Windows note: this MAGMA 1.10 build also writes `<prefix>.genes.out.txt` (extra `.txt`) instead of `<prefix>.genes.out`; `examples/magma_genebased.sh` resolves this automatically, and the same mismatch matters for PoPS below.

A tiny offline smoke-test fixture (`examples/toy_gwas_sumstats.tsv` + `examples/toy_gene_loc.txt` + `examples/toy_snp_loc.txt`, ~1,800 real 1000G-EUR SNPs across 3 gene bins spanning the real PCSK9 locus, chr1:55.4-55.6Mb hg19) ships with this Skill for sanity-checking the Steps 1-2 pipeline against a real MAGMA + PLINK reference before pointing it at real data; regenerate or rescale it with `examples/make_toy_fixture.py`.

## Open Targets L2G via GraphQL

**Goal:** Query pre-computed L2G scores for a study and locus without re-running the integrative pipeline.

**Approach:** Query the Open Targets GraphQL endpoint with a study ID and lead variant; parse per-gene scores.

### Open Targets Platform vs Genetics Portal (2024 Consolidation)

In 2024 Open Targets Genetics was merged into the Open Targets Platform GraphQL API at `api.platform.opentargets.org/api/v4/graphql`. The legacy Genetics endpoint (`api.genetics.opentargets.org/graphql`) still responds but is deprecated; new pipelines should target the Platform API. The schema also changed: the Platform exposes `credibleSet(studyLocusId)` with `l2GPredictions` (target, score, SHAP per-feature explainability), whereas the legacy schema exposed `studyLocus2GeneTable` with `yProbaModel` and per-component sub-scores.

Legacy (Genetics, deprecated):

```graphql
query L2G_legacy($studyId: String!, $variantId: String!) {
  studyLocus2GeneTable(studyId: $studyId, variantId: $variantId) {
    rows {
      gene { symbol }
      yProbaModel
      yProbaDistance
      yProbaMolecularQTL
      hasColoc
    }
  }
}
```

Modern (Platform, recommended):

```graphql
query L2G_modern($studyId: String!) {
  credibleSet(studyLocusId: $studyId) {
    l2GPredictions {
      rows {
        target { approvedSymbol }
        score
        features { name value shapValue }
        shapBaseValue
      }
    }
  }
}
```

```python
import requests
import pandas as pd

resp = requests.post('https://api.platform.opentargets.org/api/v4/graphql',
                     json={'query': '...modern L2G query...',
                           'variables': {'studyId': 'GCST006464_locus_42'}})
preds = resp.json()['data']['credibleSet']['l2GPredictions']['rows']
l2g_df = pd.json_normalize(preds).sort_values('score', ascending=False)
```

The headline `score` (Platform) corresponds to `yProbaModel` (legacy). Platform `features[].shapValue` values replace the legacy `yProba*` sub-scores and explain what drove the prediction. Genes whose SHAP is dominated by the distance feature but minimal on QTL or chromatin features are distance-only candidates; trust the integrated `score` as primary.

**L2G is populated for GWAS-type credible sets, not molecular-QTL ones.** Querying a gene's own `credibleSets` (e.g. via `target(ensemblId)`) mostly returns eqtl/pqtl/sqtl-type loci with `l2GPredictions.count == 0` -- L2G is computed against GWAS-trait credible sets, at the trait's own lead locus. Use the top-level `credibleSets(studyTypes: [gwas], ...)` query to find a GWAS-type `studyLocusId` first (verified live, 2026-09-18). `examples/opentargets_l2g_query.py` is a runnable, dependency-free (stdlib `urllib`/`json` only) template implementing exactly this two-step lookup and the modern query above; run it directly against the live, unauthenticated API.

## PoPS Polygenic Priority Score

**Goal:** Add a distance-orthogonal similarity-based prior to ranked gene candidates.

**Approach:** Run MAGMA genome-wide to produce gene Z; feed gene Z plus a curated gene-feature matrix (pathway membership + co-expression + PPI) to PoPS ridge (L2-penalized) regression. Per-gene priority scores are produced; per-locus relative ranking is informative.

```bash
# PoPS requires the gene-feature matrix and MAGMA gene Z output
# Download features and gene_annot from FinucaneLab/pops releases

python pops.py \
    --gene_annot_path gene_annot.txt \
    --feature_mat_prefix PoPS_features_full \
    --control_features_path control.features \
    --num_feature_chunks 10 \
    --magma_prefix gene_step \
    --out_prefix pops_out

# Output pops_out.preds: per-gene priority score
# Output pops_out.coefs: per-feature ridge coefficients (interpretation)
```

PoPS is biology-agnostic; the feature matrix encodes biology. Bias in the features (e.g. cancer-pathway-heavy gene sets for a non-cancer trait) propagates to the output; verify feature coverage matches the trait.

**Windows: rename MAGMA's output before running PoPS.** MAGMA 1.10 on Windows writes `<magma_prefix>.genes.out.txt` (extra `.txt`), but `pops.py` hard-codes `<magma_prefix>.genes.out` and raises `FileNotFoundError` if the exact name is missing -- following this section verbatim on Windows fails with no explanation (verified end-to-end: PoPS confirmed working immediately after `cp <prefix>.genes.out.txt <prefix>.genes.out`; corroborated independently in this Skill's tooling environment, `mendelian-randomization-analyst/TOOLS.md`). `examples/pops_run.py` wraps this exact PoPS invocation and performs the rename automatically (no-op on Linux/Mac, where MAGMA writes `.genes.out` directly) -- prefer it over calling `pops.py` directly on Windows.

**PoPS needs genome-wide, multi-chromosome MAGMA input to be meaningful.** With only one chromosome (or a single locus) of genes, PoPS's held-out-chromosome ridge CV has no held-out fold to validate against, `SELECTED_CV_ALPHA` saturates at its maximum, and every `PoPS_Score` collapses toward 0 -- the code path still runs to completion (exit 0, real output files), but the scores carry no signal at that scale (verified on a real 3-gene single-locus run). Do not present a locus-scale PoPS run's scores as a real result; use it only to confirm the pipeline is wired correctly, and require genome-wide MAGMA output for an actual PoPS-based effector-gene call.

## Multi-Evidence Integration: Concordance Scoring

**Goal:** Combine fine-mapping, coloc, ABC, L2G, PoPS, and distance into a per-locus per-gene concordance score; flag high-confidence candidates.

**Approach:** Per-locus, gather evidence per candidate gene from each method; score each evidence stream as pass / fail at the canonical threshold; sum the passing streams; report >= 3 passing as high-confidence.

```r
library(dplyr)

# Per-locus candidate gene table; one row per (locus, gene)
candidates <- read.table('locus_candidates.tsv', header = TRUE, sep = '\t')

# Score each evidence stream against canonical thresholds
candidates <- candidates %>%
  mutate(
    pass_finemap = pip_top_variant > 0.5 & credible_set_purity > 0.5,
    pass_coloc = coloc_pph4 >= 0.7,
    pass_distance = distance_to_tss <= 100000,
    pass_pops = pops_decile_rank == 1,
    pass_l2g = l2g_score >= 0.5,
    pass_abc = abc_score >= 0.02 | encode_re2g_score >= 0.5,
    concordance = pass_finemap + pass_coloc + pass_distance +
                  pass_pops + pass_l2g + pass_abc,
    confidence_tier = case_when(
      concordance >= 5 ~ 'near_certain',
      concordance >= 4 ~ 'strong',
      concordance >= 3 ~ 'high',
      concordance >= 2 ~ 'suggestive',
      TRUE ~ 'associational_only'))

# Report
candidates %>%
  filter(concordance >= 3) %>%
  arrange(desc(concordance), desc(l2g_score)) %>%
  select(locus, gene, concordance, confidence_tier, l2g_score, pops_decile_rank, coloc_pph4)
```

Concordance scoring is conservative; some real causal genes score 2-of-6 because not all evidence streams are available at all loci. Report the per-stream availability alongside the concordance score so readers know whether failure reflects negative evidence or absent evidence.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| L2G top gene != nearest gene | Distal regulation; L2G integrates fine-mapping + coloc + chromatin | Trust L2G; verify with ABC / ENCODE-rE2G if epigenome data exists |
| L2G top gene != PoPS top gene | Orthogonal methods; one is locus-feature-based, one is similarity-based | Both valid signals; concordance is high-confidence, discordance is investigate-further |
| MAGMA top gene different from L2G | MAGMA is window-based; L2G integrates distal regulation | MAGMA is a baseline; L2G is the more comprehensive scorer |
| Two genes both pass L2G >= 0.5 at the locus | Multi-effector locus or LD-tied co-regulated genes | Run FUSION.post_process.R (conditional/joint analysis by default) or coloc.susie for independence; functional validation needed |
| L2G high, no eQTL coloc | Driver may be a coding variant or chromatin-mediated regulation absent from eQTL panel | Check VEP for coding consequence; check pQTL panels (Sun 2018 Nature; UKB-PPP) |
| FUMA top gene != Open Targets L2G top gene | Different feature weighting and FUMA web-platform version may lag OT pipeline | Both informative; document which version of each |
| ABC predicts gene A, eQTL coloc predicts gene B | Cell-type mismatch in epigenome panel vs eQTL panel | Match cell type; if not possible, prefer the panel matched to causal tissue (per LDSC-SEG) |
| All methods agree | Strong concordance | Report as high-confidence; consider for CRISPRi validation |
| All methods fail at this locus | No causal gene resolvable from current data | Report as unresolved; flag for matched-tissue eQTL or epigenome data |

**Operational rule:** Concordance across L2G + PoPS + (coloc OR ABC) is the publication-grade triangulation standard. Single-method calls are weak; require >= 3 of 6 orthogonal evidence streams. Disagreement between L2G and PoPS at a locus often points to (a) a multi-effector locus, (b) a method-feature artefact, or (c) genuine biological complexity. Report both, with their per-locus ranks.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| MAGMA error `--gene-loc file required` | Path mismatch or wrong reference build | Provide hg19 or hg38 gene location file matching the GWAS build |
| MAGMA `--set-annot` errors `insufficient degrees of freedom` or `no input variables to analyse` | Too few genes in `--gene-results` for the 6-covariate competitive regression (see MAGMA pipeline section) | Skip Step 3 at locus scale; gene-set enrichment needs several hundred+ genes (genome-wide or many-loci input) |
| All MAGMA genes show p ~ 1 | LD reference panel mismatch (ancestry or sample size) | Use 1000G EUR g1000_eur reference for EUR GWAS; verify with `--gene-results` summary |
| Open Targets API returns empty | Study not in the OT release OR wrong study ID format | Verify study ID via OT study search; some studies require GCST prefix |
| PoPS output has all genes scoring near 0 | Feature matrix mis-aligned to gene_annot OR MAGMA gene Z scale wrong; OR MAGMA input is single-chromosome / locus-scale (ridge CV has no held-out fold, see PoPS section) | Verify column ordering in feature matrix; check MAGMA `genes.raw` parsing; confirm MAGMA input is genome-wide, multi-chromosome before trusting scores |
| PoPS raises `FileNotFoundError` on `<prefix>.genes.out` | Windows MAGMA wrote `<prefix>.genes.out.txt` (extra `.txt`); pops.py expects the exact name | `cp <prefix>.genes.out.txt <prefix>.genes.out` before running, or use `examples/pops_run.py` (handles this automatically) |
| cS2G gene allocation differs from L2G | Different aggregation strategies; cS2G heritability-calibrated, L2G classifier-trained | Both informative; cS2G is a per-SNP aggregator, L2G is per-(locus, gene) |
| ABC predicts a passenger gene | Wrong cell-type Hi-C or H3K27ac in ABC input | Verify cell-type-matched epigenome; cross-reference atac-seq/enhancer-gene-linking |
| FUMA SNP2GENE job stuck | Web platform queue OR exceeded GWAS size limit | Re-submit; reduce sumstats to genome-wide-significant loci if oversize |
| L2G high at HLA region | Method is not designed for long-range LD | Exclude HLA from genome-wide V2G summaries; report separately |
| Multiple candidate genes at locus, no clear winner | Multi-effector or LD-tied co-regulation | Allow multi-gene reporting; functional validation needed |
| PoPS top gene is a passenger | Feature matrix bias (pathway over-representation) | Inspect PoPS coefficients; verify feature relevance to trait |

## Tool Install Notes

- **MAGMA**: Pre-compiled binary from cncr.nl/research/magma. Linux / Mac / Windows. Ships as `magma` CLI; needs PLINK bfile reference (e.g. 1000G g1000_eur). MAGMA ships as a pre-compiled zip, not a git repo -- the old `ctg.cncr.nl` and `vu.data.surfsara.nl` direct links now redirect to an HTML landing page, so a hard-coded `wget` would save that page AS the zip and `unzip` would fail; download the current zip manually from the CNCR page, then `unzip magma_v*.zip` (filename tracks the version).
- **FUMA**: Web platform at fuma.ctglab.nl/snp2gene. No local install. Requires user account; SNP2GENE jobs run server-side with current reference annotations.
- **Open Targets Platform (current)**: GraphQL at `api.platform.opentargets.org/api/v4/graphql`. Query L2G via `credibleSet(studyLocusId)` -> `l2GPredictions { rows { target, score, features } }`. Recommended for new pipelines.
- **Open Targets Genetics (legacy, deprecated)**: GraphQL at `api.genetics.opentargets.org/graphql`. Python `pip install gentropy` (official Open Targets) or `pip install otargenpy` (community GraphQL wrapper) OR direct GraphQL queries via `requests`. Genetics Portal was consolidated into the integrated Open Targets Platform in 2024; legacy endpoint still responds but new work should target the Platform.
- **PoPS**: `git clone https://github.com/FinucaneLab/pops`. Python; ships with feature matrix download instructions. Pre-built feature matrix at the releases page.
- **cS2G**: Pre-computed gene scores downloadable from zenodo.org/records/7754032 (cS2G_UKBB.zip / cS2G_1000GEUR.zip). No install; lookup table.
- **DEPICT**: `git clone https://github.com/perslab/depict`. Java + Python; legacy method, see Pers 2015.
- **FLAMES**: Recent (Schipper M et al 2025 Nat Genet 57:323); check the publication's GitHub for the current install path.
- **INQUISIT**: Originally for breast cancer (Fachal 2020 Nat Genet 52:56); see the paper's supplementary methods for adaptation to other traits.
- **ABC**: `git clone https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction`. Python; see atac-seq/enhancer-gene-linking for the full pipeline.
- **ENCODE-rE2G**: `git clone https://github.com/EngreitzLab/ENCODE_rE2G`. Snakemake; see atac-seq/enhancer-gene-linking.

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Was the nearest gene just taken?" | No. L2G + PoPS concordance reported per locus; ABC / ENCODE-rE2G applied when matched epigenome exists. Nearest-gene flagged as low-confidence fallback only. |
| "Was L2G + PoPS concordance tried?" | Yes. Both scores reported per (locus, gene) candidate; multi-stream concordance score (>= 3 of 6 streams) gates high-confidence calls. |
| "Was the eQTL panel tissue-relevant?" | Causal tissue identified via stratified LDSC / LDSC-SEG; TWAS run in that tissue; cross-tissue check via S-MultiXcan. |
| "Was ABC / ENCODE-rE2G applied?" | Applied when matched ATAC + H3K27ac (+ Hi-C optional) available; otherwise flagged as limitation. Cross-reference atac-seq/enhancer-gene-linking. |
| "Was CRISPRi-FlowFISH validation considered?" | Cross-referenced against Fulco 2019 K562 catalog (>3,500 pairs); Gasperini 2019 at-scale catalog (~75,000 pairs); Schraivogel 2020 TAP-seq / targeted Perturb-seq screen in K562. |
| "What MAGMA window was used?" | Window stated explicitly. FUMA 35+10 vs MAGMA-native 0+0 distinction made; sensitivity over 0+0, 35+10, 50+50 reported for high-stakes loci. |

## References

- Mountjoy E, Schmidt EM, Carmona M, Schwartzentruber J, Peat G et al 2021 Nat Genet 53:1527 (Open Targets L2G; locus-to-gene integrative scorer)
- Ghoussaini M, Mountjoy E, Carmona M, Peat G, Schmidt EM et al 2021 Nucleic Acids Res 49:D1311 (Open Targets V2G/Genetics)
- de Leeuw CA, Mooij JM, Heskes T, Posthuma D 2015 PLoS Comput Biol 11:e1004219 (MAGMA gene-based association)
- Watanabe K, Taskesen E, van Bochoven A, Posthuma D 2017 Nat Commun 8:1826 (FUMA SNP2GENE integrative annotation)
- Gazal S, Weissbrod O, Hormozdiari F, Dey KK, Nasser J et al 2022 Nat Genet 54:827 (cS2G combined SNP-to-gene)
- Weeks EM, Ulirsch JC, Cheng NY, Trippe BL, Fine RS et al 2023 Nat Genet 55:1267 (PoPS Polygenic Priority Score)
- Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT et al 2019 Nat Genet 51:1664 (ABC enhancer-gene linking; CRISPRi-FlowFISH validation)
- Nasser J, Bergman DT, Fulco CP, Guckelberger P, Doughty BR et al 2021 Nature 593:238 (ABC genome-wide application)
- Gschwind AR, Mualim KS, Karbalayghareh A, Sheth MU, Dey KK et al 2023 bioRxiv 2023.11.09.563812 (ENCODE-rE2G; preprint)
- Pers TH, Karjalainen JM, Chan Y, Westra HJ, Wood AR et al 2015 Nat Commun 6:5890 (DEPICT)
- Fachal L, Aschard H, Beesley J, Barnes DR, Allen J et al 2020 Nat Genet 52:56 (INQUISIT breast-cancer V2G)
- Finucane HK, Reshef YA, Anttila V, Slowikowski K, Gusev A et al 2018 Nat Genet 50:621 (LDSC-SEG tissue prioritisation)
- Yazar S, Alquicira-Hernandez J, Wing K, Senabouth A, Gordon MG et al 2022 Science 376:eabf3041 (OneK1K sc-eQTL)
- Stacey D, Fauman EB, Ziemek D, Sun BB, Harshfield EL et al 2019 Nucleic Acids Res 47:e3 (ProGeM)
- Sun BB, Maranville JC, Peters JE, Stacey D, Staley JR et al 2018 Nature 558:73 (pQTL panel)
- Vosa U, Claringbould A, Westra HJ, Bonder MJ, Deelen P et al 2021 Nat Genet 53:1300 (eQTLGen reference)
- Mancuso N, Freund MK, Johnson R, Shi H, Kichaev G et al 2019 Nat Genet 51:675 (FOCUS for gene-level fine-mapping; cross-reference TWAS skill)
- Schaid DJ, Chen W, Larson NB 2018 Nat Rev Genet 19:491 (fine-mapping review)
- Barbeira AN, Pividori M, Zheng J, Wheeler HE, Nicolae DL, Im HK 2019 PLoS Genet 15:e1007889 (S-MultiXcan; joint per-tissue TWAS via PC decomposition)
- Hu Y, Li M, Lu Q, Weng H, Wang J et al 2019 Nat Genet 51:568 (UTMOST; cross-tissue weight imputation TWAS)
- Gasperini M, Hill AJ, McFaline-Figueroa JL, Martin B, Kim S et al 2019 Cell 176:377 (at-scale CRISPRi-FlowFISH enhancer-gene screen, ~75,000 pairs)
- Schraivogel D, Gschwind AR, Milbank JH, Leonce DR, Jakob P et al 2020 Nat Methods 17:629 (TAP-seq / targeted Perturb-seq enhancer-gene screen in K562)
- Amemiya HM, Kundaje A, Boyle AP 2019 Sci Rep 9:9354 (ENCODE blacklist v2; HLA-region exclusion reference)

## Related Skills

- causal-genomics/fine-mapping - Variant-level credible sets feeding L2G and concordance scoring
- causal-genomics/colocalization-analysis - coloc PP.H4 as one of the six evidence streams
- causal-genomics/transcriptome-wide-association - TWAS gene-level association and FOCUS fine-mapping
- causal-genomics/mendelian-randomization - cis-eQTL MR as orthogonal causal evidence
- causal-genomics/mediation-analysis - Downstream gene-mediated trait effects
- causal-genomics/proteome-mr-drug-target - pQTL-based effector-gene nomination for drug targets
- atac-seq/enhancer-gene-linking - ABC and ENCODE-rE2G enhancer-gene predictions feeding distal-regulation evidence
- atac-seq/atac-peak-calling - Generating enhancer candidates from matched-tissue ATAC
- atac-seq/deep-learning-atac - chromBPNet variant-effect prediction as 7th evidence stream for non-coding effector inference
- gene-regulatory-networks/coexpression-networks - Gene co-expression features feeding PoPS
- gene-regulatory-networks/scenic-regulons - TF-target regulons as supporting evidence at regulatory loci
- pathway-analysis/go-enrichment - Pathway context for effector-gene candidates
- population-genetics/association-testing - Upstream GWAS summary-statistic generation
- variant-calling/variant-annotation - Coding-consequence annotation for effector-gene variants
- workflows/gwas-pipeline - End-to-end GWAS pipeline producing effector-gene-prioritisation input
