# Per-Method Failure Modes: Effector-Gene Prioritization

Moved out of `SKILL.md` (bio-causal-genomics-effector-gene-prioritization) because these
per-method trigger/mechanism/symptom/fix write-ups are needed only when diagnosing a specific
discordance or pitfall -- consult this file when `SKILL.md`'s intro or Decision Tree points
you here.

## Nearest-gene assumption fails (most common pitfall)

**Trigger:** Assigning the GWAS lead variant to the closest gene without checking long-range regulation.

**Mechanism:** Approximately 30-50% of well-fine-mapped GWAS variants regulate a gene that is NOT the nearest TSS (Mountjoy 2021; Fulco 2019 Nat Genet 51:1664). Distal enhancer-promoter contacts span 50 kb to > 1 Mb; LD around the lead variant often spans only kilobases, so the credible-set centroid may sit closer to a passenger gene than to the true target.

**Symptom:** Distance-based prioritisation names the nearest gene; subsequent eQTL coloc, ABC, and ENCODE-rE2G converge on a different gene at the same locus. Functional validation (CRISPRi at the variant) confirms the distal gene.

**Fix:** Use L2G (which includes distance but does not let it dominate), PoPS (which is distance-orthogonal by construction), and ABC / ENCODE-rE2G when matched epigenome data are available. Report all candidate genes at the locus with their evidence-stream contributions; do not collapse to the nearest by default.

## eQTL tissue mis-specification

**Trigger:** Using a single-tissue eQTL panel (e.g. whole blood) when the causal tissue is different (e.g. liver for lipid traits, hypothalamus for energy balance).

**Mechanism:** Cis-eQTL effect sizes are tissue-specific; eQTLs in the wrong tissue still tag the GWAS signal via LD and produce spurious colocalisations or TWAS hits. The right gene at the wrong tissue is statistically detectable but biologically uninterpretable.

**Symptom:** Strong colocalisation in a tissue biologically irrelevant to the trait; null in the expected tissue. LDSC-SEG / CELLEX / EWCE prioritisation on the GWAS sumstats independently disagrees with the eQTL tissue.

**Fix:** Run multi-tissue eQTL coloc (e.g. all GTEx tissues via S-MultiXcan + per-tissue coloc) and prioritise the tissue identified by LDSC-SEG (Finucane 2018 Nat Genet 50:621) or CELLEX. For cell-type-specific traits, move to sc-eQTL panels (OneK1K, Yazar 2022 Science 376:eabf3041). S-MultiXcan (Barbeira 2019 PLoS Genet 15:e1007889) jointly tests per-tissue z-scores via PC-decomposition of LD-induced covariance and is preferred for standard GTEx v8 panels (pre-computed weights available); UTMOST (Hu 2019 Nat Genet 51:568) imputes cross-tissue expression weights before testing and is preferred when retraining cross-tissue weights for a custom panel.

## MAGMA gene-window choice

**Trigger:** Default `--gene-loc` with 0kb upstream / 0kb downstream window assigns all SNPs only within annotated gene bodies.

**Mechanism:** A wide window (e.g. 35kb upstream + 10kb downstream) captures more regulatory SNPs per gene but assigns each tag-SNP to multiple genes simultaneously, diluting per-gene signal and inflating false positives at gene-dense regions. A narrow window misses regulatory SNPs outside the gene body and loses true positives at intergenic enhancers. The three window conventions in circulation are not interchangeable: MAGMA-native default is 0+0 (no expansion); the MAGMA paper recommended a 50+50 sensitivity check; FUMA SNP2GENE uses 35+10 (35 kb upstream + 10 kb downstream) which is FUMA's convention, NOT MAGMA's default. For brain traits, 50+50 captures distal cis-eQTL signal; for cardiometabolic traits a tighter 10+10 is more conservative. State explicitly which window was used in methods reporting.

**Symptom:** Many genes per locus flagged at p < 0.05/22k with the wide window; few genes at all flagged with the narrow window; top genes change substantially across window choices.

**Fix:** Use a sensible default (35kb upstream + 10kb downstream is the FUMA recommendation; 0+0 is MAGMA-native; 50+50 is the MAGMA-paper sensitivity window). Always pair MAGMA with eQTL-based mapping (S-PrediXcan, coloc) for distal-regulatory signal; MAGMA alone is the lightweight baseline, not the full answer.

**Wide-window 1Mb warning:** Going to 100+100 kb or 1 Mb assigns one SNP to 8-12 genes simultaneously at gene-dense loci (e.g. MHC, chr19q13, chr17q21), diluting power and creating interpretation ambiguity. Avoid 1Mb windows; if distal regulation is suspected supplement with ABC / ENCODE-rE2G enhancer-gene linkage (cross-reference atac-seq/enhancer-gene-linking) rather than widening the MAGMA window.

## Coloc fails when the locus has multiple causal variants

**Trigger:** PP.H4 < threshold despite biological evidence that the gene is causal.

**Mechanism:** coloc.abf's single-causal-variant assumption forces posterior mass to PP.H3 (distinct causal variants) when 2+ independent signals in moderate LD drive both traits. The result is a false-negative coloc call at a true effector-gene locus.

**Symptom:** Visual LocusZoom overlap is convincing but PP.H4 stays in 0.3-0.6; coloc.susie or eCAVIAR reveals multiple credible sets and a per-credible-set PP.H4 > 0.7.

**Fix:** Run coloc.susie (not coloc.abf) at gene-dense / signal-rich loci. Cross-reference causal-genomics/colocalization-analysis; do not rely on coloc.abf as the sole coloc evidence stream when allelic heterogeneity is plausible.

## PoPS vs L2G discordance

**Trigger:** PoPS top-ranked gene at locus disagrees with L2G top-ranked gene.

**Mechanism:** PoPS uses similarity-based features (pathway membership, co-expression, PPI), L2G uses per-locus features (distance, fine-mapping, coloc, chromatin). They are orthogonal by construction; disagreement is informative, not a failure.

**Symptom:** Same locus, different top gene under each method.

**Fix:** Use BOTH and treat concordance (top gene matches across L2G and PoPS) as the strongest single-locus signal short of CRISPR validation. Concordance between PoPS and locus-based methods markedly increases positive predictive value over either method alone (Weeks 2023 Nat Genet 55:1267). Report both ranks; flag concordance.

## Pleiotropic locus / multiple causal genes per locus

**Trigger:** Two or more genes at a single GWAS locus are each independently causal (different SNPs or different mechanisms).

**Mechanism:** Standard V2G frameworks assume one causal gene per locus. Real biology violates this: an estimated 5-10% of GWAS loci have multiple causal genes (a working convention; CRISPRi-FlowFISH catalogs document multi-gene loci).

**Symptom:** Two genes at the locus both pass conditional independence checks (FUSION.post_process.R conditional/joint analysis, GCTA-COJO); both show strong eQTL coloc; both have CRISPRi support.

**Fix:** Allow multi-gene reporting. Each candidate gene needs its own credible variant set (SuSiE / coloc.susie). Report the locus as multi-effector; consider CRISPRi-FlowFISH or MPRA for ground-truth resolution. Do not force a single-gene assignment. Existing CRISPRi enhancer-gene perturbation catalogs for cross-checking computational predictions: Fulco 2019 Nat Genet 51:1664 (>3,500 CRISPRi-FlowFISH enhancer-gene connections for 30 genes in K562); Gasperini 2019 Cell 176:377 (~75,000 pairs at-scale); Schraivogel 2020 Nat Methods 17:629 (TAP-seq / targeted Perturb-seq enhancer-gene screen in K562). Cite the specific catalog when reporting "validated against CRISPRi" rather than the generic term.
