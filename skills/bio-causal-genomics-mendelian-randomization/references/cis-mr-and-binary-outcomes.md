# Drug-target / cis-MR framework and binary outcomes

Read for cis-pQTL/cis-eQTL drug-target MR, binary or case-only outcomes, and collider bias in stratified designs. Verbatim from SKILL.md.

## Drug-Target / cis-MR Framework

cis-MR restricts instruments to the cis-regulatory window of the gene encoding the protein/transcript exposure (Schmidt 2020 Nat Commun 11:3255), relaxing the exclusion-restriction assumption because the protein product directly mediates the SNP's effect on the outcome. Operational core: extract cis-pQTL/cis-eQTL within +/-500 kb of the gene; clump at r2 < 0.1 (looser than polygenic MR to retain power within a narrow window); require colocalization PP.H4 >= 0.7; flag protein-altering variants (PAV) which can break SomaScan/Olink aptamer/antibody binding rather than reflect biology.

Full drug-target cis-MR workflow including UKB-PPP / deCODE / Fenland pQTL panels, PAV flagging, Olink vs SomaScan replication (~15-30% cross-platform disagreement), and the operational claim ladder lives in causal-genomics/proteome-mr-drug-target. Use that skill for any drug-target nomination.

### Binary outcomes and non-collapsibility

MR with logistic-GWAS sumstats returns per-allele log-OR on the **population-averaged** (marginal) scale, NOT the conditional log-OR (Burgess 2017 Stat Methods Med Res 26:2333). For rare disease (prevalence < 10%), OR ~= RR ~= HR and the distinction is harmless. For common disease, OR diverges from RR/HR and the MR estimate cannot be back-converted to a conditional effect without strong assumptions; report as "per 1-SD increase in genetically-predicted X, OR for Y = ..." rather than implying an individual-level intervention effect.

Collider bias when conditioning on a collider variable (Coscia 2022 Eur J Epidemiol 37:671 formalizes this for stratified MR): case-only or disease-progression designs condition the sample on disease status, opening a collider path between any cause of disease and any cause of progression. MR within affected subsets without explicit adjustment for selection probability is fragile; weight by inverse probability of selection or restrict claims to the unconditioned population.
