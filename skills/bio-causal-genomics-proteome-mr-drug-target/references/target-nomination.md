# Drug Repurposing and Target Nomination

Moved verbatim from SKILL.md. Read when nominating a target for repurposing or cross-checking against Open Targets, CTD and DGIdb.

## Drug Repurposing and Target Nomination

The Open Targets Drug platform (Ochoa 2021 Nucleic Acids Res 49:D1302) integrates approved-drug-target relationships, cis-pQTL/cis-eQTL MR, and locus-to-gene (L2G) scores (Mountjoy 2021). A target is nominated for repurposing when:

- L2G score at the GWAS lead points to the target gene
- Cis-MR estimate concordant with disease direction (lower protein -> lower disease for inhibitor candidate)
- Coloc PP.H4 >= 0.7 between target's cis-pQTL and disease GWAS
- A licensed drug exists that modulates the target
- The on-target adverse-effect pheWAS is acceptable

The PCSK9 monoclonal-antibody story is the canonical positive example; the CETP-inhibitor story is a canonical cautionary tale (cis-MR underestimated trial result due to off-target effects).

Cross-validate target nominations against the Comparative Toxicogenomics Database (CTD; ctdbase.org) and the Drug-Gene Interaction Database (DGIdb; dgidb.org) for established drug-target evidence and tractability annotations. STROBE-MR (Skrivankova 2021 JAMA 326:1614) provides the 20-item checklist for drug-target MR reporting; see causal-genomics/pleiotropy-detection for the full table.

