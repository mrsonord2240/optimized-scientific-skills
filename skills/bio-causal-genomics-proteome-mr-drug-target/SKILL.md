---
name: bio-causal-genomics-proteome-mr-drug-target
description: Runs cis-pQTL Mendelian randomization for drug-target validation using UKB-PPP (Olink), deCODE (SomaScan), Fenland, INTERVAL, ARIC, and FinnGen-PPP proteomes plus colocalization triangulation, phenome-wide on-target adverse-effect scans, cross-platform Olink/SomaScan replication, and PAV (protein-altering variant) sensitivity. Use when nominating or de-risking a drug target from plasma-proteome GWAS, mimicking pharmacological inhibition via cis-pQTL instruments, separating shared-causal from LD-confounded signal under the Schmidt 2020 cis-MR framework, screening on-target adverse phenotypes pheWAS-style, or producing publication-grade STROBE-MR plus PP.H4 evidence for a target gene.
tool_type: mixed
primary_tool: TwoSampleMR
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: TwoSampleMR 0.5.11+, MendelianRandomization 0.10+, MR-PRESSO 1.0+, coloc 5.2.3+, susieR 0.12.35+, ieugwasr 1.0+, plink2 2.00a5+, R 4.4+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `plink2 --version`; VEP `vep --help`

If code throws OAuth or rate-limit errors from OpenGWAS, or a missing `dataset$N` from coloc, introspect the installed API and adapt the example rather than retrying. UKB-PPP, deCODE, and Fenland summary statistics changed file layouts between 2023 and 2025; verify column headers before passing into `format_data()`.

# Proteome-Wide Drug-Target Mendelian Randomization

**"Does genetically lowering plasma protein X cause a change in disease Y, mimicking a drug?"** -> Use cis-pQTLs in the gene window for protein X as instruments under the Schmidt 2020 framework (Nat Commun 11:3255), restrict the exclusion-restriction violation to the geometric neighbourhood of the encoding gene, triangulate with colocalization (3-tier PP.H4 ladder below) and cross-platform replication (Olink vs SomaScan), and flag protein-altering-variant (PAV) confounding. A single significant cis-MR estimate is necessary but not sufficient for a drug-target claim; the operational bar is MR + coloc + cross-platform agreement + PAV-excluded sensitivity.

### PP.H4 Three-Tier Threshold Ladder

| Tier | PP.H4 | Use case |
|------|-------|----------|
| Suggestive | >= 0.7 | Open Targets Genetics (Mountjoy 2021 Nat Genet 53:1527) / exploratory; consistent with shared-causal |
| Standard publication | >= 0.8 | Wallace 2020 PLoS Genet 16:e1008720; most peer-reviewed pubs |
| Industry / clinical | >= 0.95 | Drug-claim grade; pharma internal target-validation standard (clinical-pharmacology bar) |

Operational rule: drug-target nomination requires PP.H4 >= 0.8 minimum; industry-grade clinical claim requires PP.H4 >= 0.95 plus the full triangulation panel.

- R (canonical): `TwoSampleMR::mr()` orchestrates the cis-IVW + Egger + median + Wald-ratio panel
- R (correlated cis-pQTLs in a window): `MendelianRandomization::mr_input(..., correlation = ld_matrix)` then `MendelianRandomization::mr_ivw(mr_obj, model = 'default')` (namespace it; see `references/correlated-instruments.md`)
- R (triangulation): `coloc::coloc.abf()` or `coloc::coloc.susie()` on the same cis-window
- pheWAS: `ieugwasr::associations()` against the OpenGWAS catalogue, looped over outcomes
- VEP CLI: annotate every cis-pQTL with `vep --species homo_sapiens --canonical --check_existing` for PAV flagging

## Reference Files

| File | Read when |
|------|-----------|
| `references/pqtl-datasets.md` | Choosing or citing a pQTL source (UKB-PPP, deCODE, Fenland, INTERVAL, ARIC, FinnGen-PPP), checking platform coverage, pinning versions |
| `references/failure-modes.md` | Platform discordance and the Olink panel pre-check, neighbour-gene LD, reverse causation (Steiger), PAV confound and the concordance rule, trans-pQTLs, UKB sample overlap |
| `references/correlated-instruments.md` | Several cis-pQTLs in moderate LD: `mr_input(correlation=)`, `MendelianRandomization::mr_ivw` namespace and sign alignment, robust/penalized IVW, Patel 2023 |
| `references/phewas.md` | Phenome-wide on-target adverse-effect scan (curated-endpoint loop, Bonferroni denominator) |
| `references/target-nomination.md` | Repurposing nomination criteria, Open Targets, CTD, DGIdb |

## Cis-MR Methodological Taxonomy

| Method | Cis-window assumption | Min cis-pQTLs | Strength | Fails when |
|--------|------------------------|----------------|----------|------------|
| Single cis-pQTL Wald ratio | Single sentinel SNP within +/-500 kb | 1 | Simplest, transparent point estimate `beta_Y/beta_X` and ratio SE | Confounded by LD-linked eQTL/pQTL of neighbour gene; no heterogeneity test |
| Cis-IVW (clumped r2 < 0.1) | Multiple weakly-correlated cis-pQTLs | 2 | Pools information, increases precision (Schmidt 2020) | r2 between pQTLs > 0.1 inflates SE under independence assumption |
| Cis-IVW with LD matrix (`MRInput` correlation) | Cis-pQTLs in moderate LD; supply LD matrix | 2 | Correct SE under correlated instruments (Burgess, Zuber, Valdes-Marquez, Sun, Hopewell 2017 Genet Epidemiol 41:714-725) | LD matrix mismatched to summary-stat ancestry |
| Cis-MR-Egger | Directional pleiotropy across cis-pQTLs | 3+ (>=10 for power) | Sensitivity for in-window directional pleiotropy | Underpowered <10 cis-pQTLs; NOME violation `I^2_GX < 0.9` |
| Cis-weighted-median | Up to 50% invalid cis-pQTLs | 3+ | Robust to a minority of bad instruments | >50% invalid cis-pQTLs |
| MR-PRESSO in cis-window | Outlier cis-pQTLs from LD-confounded neighbours | 4+ | Removes neighbour-eQTL-tagged cis-pQTLs; distortion test (Verbanck 2018) | <4 instruments; underpowered global test |
| Conditional cis-MR with weak genetic factors (Patel 2023 Biometrics 79:3458) | Factor-analysis dimension reduction of correlated cis-pQTLs + weak-factor-robust conditional inference | 2+ | Methodologically modern alternative when instruments are highly correlated | Newer; benchmarks evolving; separate implementation from mr_ivw |
| Generalized cis-IVW with correlated SNPs | Joint multivariable cis-window | 2+ | Sound under any LD provided matrix supplied | Numerical instability when r2 ~ 1 (collinear) |
| coloc.susie + Wald-ratio per CS | Multiple independent cis-signals (allelic heterogeneity) | 1+ per credible set | Per-signal MR + per-signal PP.H4 | Requires ancestry-matched LD; spurious CS under mismatch |

Verify against Burgess 2023 *Wellcome Open Res* "Guidelines for performing Mendelian randomization" (v3+) and the Open Targets Genetics drug-target pipeline (Mountjoy 2021) before pinning a primary method for a publication.

## Decision Tree by Scenario

| Scenario | Primary method | Triangulation | Why |
|----------|----------------|----------------|-----|
| Single drug target, single sentinel cis-pQTL, single outcome | Wald ratio | coloc.abf PP.H4 + cross-platform replication + PAV flag | Minimum publishable cis-MR; simplest and most transparent |
| Single drug target, multiple independent cis-pQTLs, single outcome | Cis-IVW with in-window LD matrix | coloc.susie per credible set + cross-platform | Pools signal; correctly handles within-window LD (`references/correlated-instruments.md`) |
| Drug target with secondary independent cis-signal (allelic heterogeneity) | coloc.susie + per-CS Wald ratio | Compare effect direction across CSs | Each independent signal is its own instrument; report each |
| Phenome-wide MR on a single target | Loop Wald ratio or cis-IVW across hundreds of outcomes | Bonferroni over outcomes; coloc PP.H4 on top hits | On-target adverse-effect discovery (e.g. PCSK9 -> T2D); see `references/phewas.md` |
| On-target adverse-effect scan (already-marketed drug) | pheWAS cis-MR vs all FinnGen / OpenGWAS phenotypes | Bonferroni + coloc + clinical-event registry | Re-derives known and novel on-target effects |
| Cross-platform replication required for clinical claim | Run independently on UKB-PPP (Olink) AND deCODE (SomaScan) | Direction agreement + magnitude within 2x | Single platform never sufficient for therapeutic decision (`references/failure-modes.md`) |
| Trans-pQTL "wants" to be an instrument | Refuse | -- | Trans = horizontal pleiotropy by definition; use only as confirmatory (`references/failure-modes.md`) |
| Target gene not on Olink panel | deCODE / Fenland SomaScan only | Cross-replicate across two SomaScan cohorts | Olink panel is gated; do not infer "untested" as null (`references/pqtl-datasets.md`, panel pre-check in `references/failure-modes.md`) |
| Target in cis-region with strong neighbour eQTL | coloc.susie + coloc to non-target gene's pQTL/eQTL | Drop cis-pQTLs that coloc with non-target | Mandatory: must rule out neighbour-gene mediation (`references/failure-modes.md`) |
| Sample overlap (UKB-PPP exposure + UKB phenotype outcome) | MR-RAPS with one-sample-equivalent correction OR independent outcome (FinnGen, BBJ) | Repeat in non-overlapping cohort | Pretending overlap is two-sample inflates effect estimates (`references/failure-modes.md`) |

## Per-Method Failure Modes

Six failure modes, each with trigger, mechanism, symptom and fix in `references/failure-modes.md`: Olink vs SomaScan platform discordance (with the Olink panel coverage pre-check), LD-based pleiotropy in the cis-window, reverse causation (Steiger), PAV confound (with the PAV-excluded concordance rule), trans-pQTL pleiotropy, and sample overlap when both ends are UKB.

## Triangulation Requirement (Operational Postdoc Rule)

A drug-target causal claim that survives peer review and informs pharmacology requires MULTIPLE concordant streams, not a single significant cis-MR p-value:

1. **Cis-MR estimate** with `P < 0.05 / N_proteins` (Bonferroni proteome-wide ~ 1.7e-5 for 2923 Olink proteins) OR `P < 0.05 / N_outcomes` (Bonferroni pheWAS-wide) -- depending on the testing regime
2. **Colocalization PP.H4 >= 0.8** between cis-pQTL and outcome GWAS at the gene locus -- standard publication tier; >= 0.95 for industry-grade claim (cross-reference causal-genomics/colocalization-analysis)
3. **Cross-platform replication** -- significant on both Olink (UKB-PPP) AND SomaScan (deCODE or Fenland); direction agreement is mandatory, magnitude within 2x is acceptable
4. **Cross-cohort replication** (ideal but not strictly required) -- UKB-PPP -> deCODE -> FinnGen-PPP step-up
5. **PAV-excluded sensitivity** -- the cis-MR survives when missense / nonsense / splice / aptamer-binding-region SNPs are dropped
6. **No neighbour-gene coloc** -- no cis-pQTL in the instrument set coloc-shares a causal variant with any non-target gene's eQTL or pQTL within +/-500 kb
7. **Open Targets L2G concordance** -- Open Targets Platform locus-to-gene score for the target gene >= 0.5 at the disease GWAS lead (cross-reference causal-genomics/effector-gene-prioritization)

Operational claim ladder: cis-MR significant alone = exploratory; +coloc PP.H4 >= 0.7 = consistent with shared-causal; +cross-platform = consistent across detection chemistries; +PAV-excluded + neighbour-clear = publication-grade target nomination; +cohort replication + Open Targets L2G >= 0.5 = clinical-pharmacology-grade. Clinical-pharmacology claims require ALL 6 original criteria PLUS L2G concordance.

Report (STROBE-MR): cis-MR estimate + 95% CI, PP.H4, PAV-excluded estimate, platform agreement, neighbour-gene coloc results, sample-overlap statement, and the claim-ladder rung reached.

## Phenome-Wide Drug-Target MR

Hold the cis-pQTL instrument set fixed, loop outcomes over a curated endpoint list (FinnGen DF12, Open Targets trait map or a phecode hierarchy), Bonferroni over the number of endpoints. Rationale in `references/phewas.md`; run `Rscript scripts/phewas_curated_endpoints.R <pqtl.tsv> <endpoints.tsv> <out.tsv>` (needs an OpenGWAS token).

## Cis-MR Standard Workflow

**Goal:** Produce a defensible cis-MR estimate with triangulation, given a target gene and an outcome.

**Approach:** Extract cis-pQTLs in +/-500 kb of the gene -> compute F per instrument from exposure -> clump within window at r2 < 0.1 -> harmonise with outcome -> run TwoSampleMR -> run coloc.abf on the same window -> PAV-annotate via VEP -> report panel.

The complete workflow is `examples/cis_pqtl_mr.R` (cis-window extraction, F filter, PAV flag, clumping, harmonisation, MR panel, correlated IVW, coloc.abf, PAV-excluded rerun, decision summary). Edit the target block at its top and the `data/` paths, then run `Rscript examples/cis_pqtl_mr.R`.

## Cis-IVW with Correlated Instruments

When cis-pQTLs are correlated (r2 0.1 to 0.7), supply an ancestry-matched LD matrix via `mr_input(..., correlation = ld)` and call `MendelianRandomization::mr_ivw()` (always namespaced; `TwoSampleMR::mr_ivw` masks it). Pre-prune at r2 < 0.95. Decision rule, sign alignment and the robust/penalized and Patel 2023 alternatives are in `references/correlated-instruments.md`.

## PAV Annotation via VEP

**Goal:** Tag every cis-pQTL with its most severe coding consequence to enable a PAV-excluded sensitivity panel.

**Approach:** Format SNPs as VEP input, run VEP, parse the `Consequence` column.

```bash
echo -e "chr1\t55039548\t.\tG\tT" > cis_pqtls.vcf
vep --species homo_sapiens --assembly GRCh38 --canonical --check_existing \
    --input_file cis_pqtls.vcf --output_file pqtl_vep.tsv --tab --force_overwrite
```

PAV consequences to flag (drop in sensitivity analysis): `missense_variant`, `stop_gained`, `stop_lost`, `frameshift_variant`, `splice_acceptor_variant`, `splice_donor_variant`, `start_lost`, `protein_altering_variant`. For aptamer panels, additionally consider `synonymous_variant` within the SOMAmer-binding region (rare but documented).

## Cis-Window Width: 500 kb vs 1 Mb Decision

- **Default ±500 kb** (Schmidt 2020 Nat Commun 11:3255 standard) -- balances cis-specificity against power; matches Open Targets Genetics defaults
- **Widen to ±1 Mb** when: (a) target gene has a documented distal regulatory element in ENCODE-rE2G or ABC enhancer-gene maps; (b) < 2 genome-wide-significant pQTLs are present in ±500 kb; (c) the target gene is unusually large (gene body > 500 kb itself, e.g. DMD, RBFOX1)
- **Narrow to ±250 kb** when high-density cis-eQTL background causes multi-gene pleiotropy concerns (e.g. HLA region; gene-dense pericentromeric loci)

Pre-specify the window in the methods section; widen to ±1 Mb only with a documented distal-regulatory rationale, since window-width sensitivity is a recognized peer-review pushback.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Bonferroni for cis-MR pheWAS | Standard | `P < 0.05 / N_outcomes` for on-target adverse-effect scan |
| Coloc PP.H4 | Wallace 2020 PLoS Genet 16:e1008720 | Three-tier ladder at the top of this file |
| Bonferroni for proteome-wide cis-MR | Standard | `P < 0.05 / N_proteins` ~ 1.7e-5 for 2923 Olink proteins; ~1e-5 for 4907 SomaScan |
| Cis-pQTL F >= 10 | Staiger & Stock 1997; Burgess 2011 | Weak-instrument floor |
| r2 < 0.1 clumping in cis-window | Schmidt 2020 | Reduces LD-based pleiotropy while retaining power; looser than the r2 < 0.001 polygenic-MR convention because the window itself is the pruning mechanism |
| N >= 2 pQTL datasets in agreement | Best-practice (UKB-PPP + deCODE) | Cross-platform replication mandatory for clinical claim |
| PAV-excluded sensitivity | Sun 2023 supplementary | Required for clinical claim; Olink/SomaScan vulnerable to PAV artifact |
| Neighbour-gene coloc PP.H4 < 0.5 | Operational | Cis-pQTL must NOT coloc with non-target gene; drop if it does |
| Steiger p > 0.05 (correct direction) | Hemani 2017 PLoS Genet 13:e1007081 | Directionality check; subject to Lutz 2022 caveat |
| Sample-overlap correction if both ends UKB | Burgess 2016 Genet Epidemiol 40:597 | One-sample-equivalent bias correction |

## Reconciliation: When Evidence Streams Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Cis-MR significant on Olink, null on SomaScan | Platform-specific epitope/aptamer artifact | PAV-annotate; flag protein as platform-discordant; do not claim |
| Cis-MR sig, coloc PP.H4 < 0.5 | LD-confounded signal; not shared causal | Downgrade; cis-MR alone insufficient; do not claim drug target |
| Coloc PP.H4 high, cis-MR null | Underpowered cis-MR (few cis-pQTLs, weak F) OR shared-causal for non-causal protein | Inspect cis-pQTL strength; consider increasing window to 1 Mb |
| Cis-MR sig with all pQTLs, null after PAV exclusion | PAV artifact dominating instrument | Report PAV-excluded as primary; original as supplementary |
| Cis-MR sig in UKB-PPP -> UKB outcome, null in FinnGen outcome | Sample-overlap one-sample bias | Treat FinnGen as truth; report UKB-on-UKB as biased upward |
| Cis-pQTL coloc with non-target gene's eQTL (PP.H4 >= 0.5) | Neighbour-gene mediation | Drop this cis-pQTL; re-run cis-MR with clean instruments |
| Two independent cis-pQTLs give opposite direction effects | Allelic heterogeneity with distinct biology | Run coloc.susie per credible set; report each signal separately |
| Wald ratio at sentinel SNP differs from cis-IVW | One outlier cis-pQTL dominates IVW | Run MR-PRESSO; check Egger intercept |

**Operational rule for publication:** the claim ladder under Triangulation Requirement applies unchanged; Wald ratio replaces cis-IVW when N=1 SNP. Any single missing leg downgrades the claim to "consistent with" rather than "evidence for."

## Drug Repurposing and Target Nomination

Nominate for repurposing only with L2G pointing at the target, concordant cis-MR direction, PP.H4 >= 0.7, a licensed modulator and an acceptable on-target pheWAS. Full criteria and cross-checks (CTD, DGIdb, STROBE-MR) in `references/target-nomination.md`.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `harmonise_data` drops most SNPs | EAF columns missing or palindromic at MAF~0.5 | Provide EAF; use `action = 2` (default) or `action = 3` for strictest |
| F-statistic from outcome | Computed `beta.outcome / se.outcome` | F must come from exposure (Burgess 2011) |
| OpenGWAS OAuth failure | Token expired (OpenGWAS auth tightened 2024) | Use local plink + 1KG bfile via `ieugwasr::ld_clump(..., bfile=...)` |
| `mr_ivw()` throws `unused arguments (model=, correl=)` | `TwoSampleMR::mr_ivw` masks `MendelianRandomization::mr_ivw` | Call `MendelianRandomization::mr_ivw()` explicitly (`references/correlated-instruments.md`) |
| Cis-IVW null but Wald ratio at lead significant | Inclusion of weak / outlier cis-pQTLs | Tighten clumping; run MR-PRESSO outlier test |
| Coloc PP.H3 dominant | Multiple causal in moderate LD | Switch to coloc.susie with ancestry-matched LD |

## Tool Installation Notes

```r
install.packages(c('remotes', 'MendelianRandomization', 'coloc', 'susieR', 'dplyr'))
remotes::install_github('MRCIEU/TwoSampleMR')
remotes::install_github('MRCIEU/ieugwasr')
remotes::install_github('MRCIEU/genetics.binaRies')   # bundles plink binary
remotes::install_github('rondolab/MR-PRESSO')
```

```bash
# Ensembl VEP for PAV annotation
conda install -c bioconda ensembl-vep
vep_install -a cf -s homo_sapiens -y GRCh38 -c $HOME/.vep

# 1000 Genomes EUR plink reference for local clumping / LD matrix (bfile prefix '1kg_EUR/EUR' in `examples/cis_pqtl_mr.R` and `references/correlated-instruments.md`).
# Same panel the OpenGWAS API uses, per the ieugwasr local-LD vignette (link live 2026-09-21, ~1.5 GB).
# Use this when OpenGWAS needs a token or is unreachable.
curl -O http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
tar -xzf 1kg.v3.tgz             # yields EUR.bed/.bim/.fam (and AFR, AMR, EAS, SAS)
mkdir -p 1kg_EUR && mv EUR.* 1kg_EUR/
```

Local `ld_clump()` / `ld_matrix()` need only that bfile plus `plink_bin`; no OpenGWAS call is made when `bfile=` is given. Vignette: https://mrcieu.github.io/ieugwasr/articles/local_ld.html.

pQTL data acquisition: UKB-PPP via the UK Biobank pre-published portal (ukb-ppp.gwas.eu, per protein by UniProt ID); deCODE via https://www.decode.com/summarydata/ with a data-use agreement; Fenland via the EBI GWAS catalog and Pietzner 2021 supplementary; FinnGen-PPP via the FinnGen DF12 release portal; OpenGWAS hosts many pre-formatted pQTL studies but always verify the upstream reference and download date.

## References

- Schmidt AF et al 2020 Nat Commun 11:3255 (cis-MR drug-target framework)
- Sun BB et al 2023 Nature 622:329 (UKB-PPP Olink Explore)
- Ferkingstad E et al 2021 Nat Genet 53:1712 (deCODE SomaScan)
- Pietzner M et al 2021 Science 374:eabj1541 (Fenland SomaScan)
- Sun BB et al 2018 Nature 558:73 (INTERVAL SomaScan)
- Zhang J et al 2022 Nat Genet 54:593 (ARIC multi-ancestry pQTL)
- Emilsson V et al 2018 Science 361:769 (AGES-Reykjavik)
- Eldjarn G et al 2023 Nature 622:348 (Olink vs SomaScan cross-platform comparison)
- Schmidt AF et al 2017 Lancet Diabetes Endocrinol 5:97 (PCSK9 -> T2D pheWAS exemplar)
- Burgess S, Zuber V, Valdes-Marquez E, Sun BB, Hopewell JC 2017 Genet Epidemiol 41:714-725 (correlated-IV IVW)
- Burgess S et al 2016 Genet Epidemiol 40:597 (sample-overlap correction)
- Mountjoy E et al 2021 Nat Genet 53:1527 (Open Targets Genetics L2G + coloc)
- Ochoa D et al 2021 Nucleic Acids Res 49:D1302 (Open Targets Drug platform)
- Lutz SM et al 2022 Genet Epidemiol 46:139 (Steiger caveat under unmeasured confounding)
- Verbanck M et al 2018 Nat Genet 50:693 (MR-PRESSO)
- Wallace C 2020 PLoS Genet 16:e1008720 (coloc p12 sensitivity)
- Skrivankova VW et al 2021 JAMA 326:1614 (STROBE-MR)
- Patel A, Gill D, Newcombe P, Burgess S 2023 Biometrics 79:3458-3471 (robust cis-MR with correlated instruments)
- Mounier N & Kutalik Z 2023 Genet Epidemiol 47:314 (MRlap sample-overlap + winner's-curse correction)

## Related Skills

- causal-genomics/mendelian-randomization - Parent polygenic-MR framework; cis-MR is the drug-target specialization; MRlap sample-overlap correction
- causal-genomics/colocalization-analysis - Required PP.H4 triangulation for any cis-MR drug-target claim
- causal-genomics/fine-mapping - Credible-set construction prior to coloc.susie at the cis-locus
- causal-genomics/pleiotropy-detection - MR-PRESSO / Egger diagnostics adapted to cis-window; STROBE-MR 20-item checklist
- causal-genomics/transcriptome-wide-association - eQTL-based parallel evidence for the same target
- causal-genomics/mediation-analysis - Step from cis-MR to downstream mediator pathway
- causal-genomics/effector-gene-prioritization - Open Targets L2G concordance leg of the triangulation panel
- population-genetics/association-testing - Source GWAS pipelines for pQTL discovery
- population-genetics/linkage-disequilibrium - LD-matrix construction for correlated cis-IVW
- variant-calling/variant-annotation - VEP PAV annotation for sensitivity analysis
- clinical-databases/clinvar-lookup - Pathogenic-variant context for nominated targets
