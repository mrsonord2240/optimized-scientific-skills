# Per-Method Failure Modes

Moved verbatim from SKILL.md. Read when a cis-MR result looks platform-specific, LD-confounded, reverse-causal, PAV-driven, trans-instrumented or sample-overlapped, and for the Olink panel pre-check.

## Per-Method Failure Modes

### Olink vs SomaScan platform discordance

**Trigger:** A cis-pQTL effect size or even direction differs between UKB-PPP (Olink antibody) and deCODE/Fenland (SomaScan aptamer) for the same protein.

**Mechanism:** Olink uses paired antibody proximity-extension assay (PEA) binding distinct epitopes; SomaScan uses single-aptamer SOMAmer binding a folded epitope. A missense SNP that alters one epitope produces an apparent pQTL on that platform only (a pseudo-PAVQTL / aptamer-affinity QTL, AAVQTL). Splice/isoform differences also produce platform-specific signal. A direct Olink-vs-SomaScan cross-platform comparison reported only modest cross-platform correlation, with roughly half of cis-pQTL signals not shared across platforms (72% of Olink vs 43% of SomaScan assays had a cis-pQTL) (Eldjarn G et al 2023 Nature 622:348).

**Symptom:** Cis-MR significant on one platform, null on the other; or significant on both but opposite direction.

**Fix:** Replicate every cis-MR claim on at least one alternate platform; annotate cis-pQTLs with VEP and flag missense / nonsense / splice variants in the gene; consult Olink and SomaScan documentation for the protein's antibody / aptamer binding region; perform a PAV-excluded sensitivity analysis and report both estimates. If platforms disagree irreconcilably, report the protein as platform-discordant and do NOT advance for clinical claim.

**Olink panel coverage pre-check:** Olink Explore 3072 covers ~3,000 proteins; Explore Expansion ~3,000; Explore HT ~5,400. Always confirm the target is on the panel BEFORE running cis-MR: `grep -i <GENE> olink_panel_proteins.tsv` (panel manifest from olink.com). Sun 2023 Nature supplementary Table S1 lists every UKB-PPP-covered protein explicitly; absence from Table S1 means the target was not measured in Phase 1 (regardless of biology) and the analysis must use SomaScan.

### LD-based pleiotropy in the cis-window

**Trigger:** A cis-pQTL for target gene X is also a strong eQTL or pQTL for a neighbouring gene within +/-500 kb.

**Mechanism:** The instrument's effect on outcome may be mediated by the neighbour gene's protein, not by target X. Cis-window proximity does NOT guarantee specificity.

**Symptom:** Colocalization with the non-target gene's pQTL or eQTL returns PP.H4 >= 0.5; cis-MR using only "clean" cis-pQTLs (those not coloc'd with neighbours) gives a substantially different estimate.

**Fix:** For every cis-pQTL, run colocalization against all eQTL/pQTL signals within +/-500 kb in the relevant tissue; drop cis-pQTLs that coloc (PP.H4 >= 0.5) with any non-target gene; coloc.susie is preferred when multiple credible sets exist in the window. Reports must list which cis-pQTLs were retained and the rationale.

### Reverse causation from disease state on plasma protein

**Trigger:** The outcome trait elevates the protein as a downstream consequence (e.g. CRP elevated in coronary disease patients; TNF in autoimmune disease).

**Mechanism:** Plasma proteomes measured in observational cohorts (including UKB-PPP) include people who have or will develop the outcome; the observed protein-disease association may be downstream not upstream.

**Symptom:** Observational protein-disease association is large; cis-MR estimate is much smaller, null, or in the opposite direction.

**Fix:** Apply Steiger filtering on each cis-pQTL (`steiger_filtering()`); restrict to pre-symptomatic samples where possible (pediatric or early-adult cohorts); replicate in longitudinal cohorts measuring protein years before disease onset. Note: Steiger has its own caveat under unmeasured confounding (Lutz SM et al 2022 Genet Epidemiol 46:139); cross-validate via bidirectional cis-MR.

`steiger_filtering()` needs `samplesize.exposure`/`samplesize.outcome` on `dat` (or, for a binary/log-odds outcome, `ncase.outcome`/`ncontrol.outcome`/`prevalence.outcome` instead). Neither `format_data()` nor `read_outcome_data()` add sample size by default -- including in this Skill's own `examples/cis_pqtl_mr.R` -- and without it `add_rsq()` silently produces no `rsq`/`effective_n` column, then `steiger_filtering()` crashes (`replacement has 0 rows, data has N`; checked on TwoSampleMR 0.7.9). Set it explicitly first.

```r
library(TwoSampleMR)
dat <- harmonise_data(exposure_pQTL, outcome_GWAS)
dat$samplesize.exposure <- 54219    # UKB-PPP N here; use the real exposure GWAS N
dat$samplesize.outcome <- 122733    # outcome GWAS N (e.g. CARDIoGRAMplusC4D)
dat <- steiger_filtering(dat)
dat_forward <- dat[dat$steiger_dir, ]   # drop reverse-direction SNPs
dir_test <- directionality_test(dat_forward)
```

### PAV (protein-altering-variant) confound

**Trigger:** A cis-pQTL is itself a missense, nonsense, frameshift, splice-site, or stop-gain variant in the target gene.

**Mechanism:** The variant changes the protein sequence, which can change antibody affinity (Olink) or SOMAmer affinity (SomaScan) without changing the actual protein abundance in plasma. The pQTL appears strong but reflects measurement artifact.

**Symptom:** The strongest cis-pQTL is a coding variant; cis-MR effect magnitude shrinks substantially when PAVs are excluded; the pQTL is platform-specific.

**Fix:** Annotate ALL cis-pQTLs with Ensembl VEP (`--check_existing --canonical`). Tabulate every cis-pQTL's most-severe consequence. Run two cis-MR analyses: (a) all cis-pQTLs, (b) PAV-excluded. Report both; require concordance for a publishable claim (Sun 2023 supplementary). For aptamer panels, also annotate the SOMAmer binding-region overlap if available. Some proteins (complement factors, immunoglobulins) carry many cis PAVs; the PAV-excluded panel may discard most instruments, and then no cis-MR claim can be made.

**PAV-excluded concordance rule:** Concordance between all-cis and PAV-excluded estimates requires (i) effect direction preserved, (ii) |effect| within 2x of the all-cis estimate, AND (iii) p-value still nominally significant (P < 0.05) after PAV exclusion. If ANY of the three criteria fails, report both estimates and downgrade the claim from "drug-target" to "suggestive cis association requiring orthogonal confirmation."

### Trans-pQTL pleiotropy if used as instrument

**Trigger:** Including trans-pQTLs (outside +/-500 kb of the gene) in the instrument set to gain power.

**Mechanism:** A trans-pQTL acts through some other gene's protein that then regulates target X; using it as an instrument violates the exclusion-restriction by definition (horizontal pleiotropy).

**Symptom:** Effect estimate shifts when trans-pQTLs are added; Egger intercept significant.

**Fix:** Restrict instrument set to cis-window only (Schmidt 2020). Trans-pQTLs may be confirmatory ("does the regulator gene also predict outcome?") but never primary instruments for a drug-target claim.

### Sample overlap when both ends are UKB

**Trigger:** Using UKB-PPP for the exposure (Olink protein) AND a UKB phenotype (HES, cancer registry, ICD-10) for the outcome.

**Mechanism:** The same individuals contribute to both summary statistics; weak-instrument bias is now one-sample-equivalent and points TOWARD the confounded observational estimate, not the null.

**Symptom:** UKB-on-UKB cis-MR effect is substantially larger than the same protein-disease pair tested in an independent cohort.

**Fix:** Use an independent outcome cohort whenever possible (FinnGen, Biobank Japan, MVP). If UKB-on-UKB is necessary, apply MR-RAPS with one-sample-equivalent treatment OR the Burgess 2016 (Genet Epidemiol 40:597) sample-overlap correction. Document the overlap fraction. MRlap (Mounier 2023 Genet Epidemiol 47:314) is the recommended modern tool for joint sample-overlap and winner's-curse correction; see causal-genomics/mendelian-randomization for the canonical implementation.

