# Input-dependent alternatives
LCV, LHC-MR and bidirectional MR are source-described alternatives, but supplied examples depend on study-specific GWAS inputs, preprocessing and/or pre-existing objects. They are not turnkey file-based runs in this draft. Stop and adapt a validated workflow with appropriate ancestry, LD references, sample-overlap handling and harmonization.
- [LCV](../references/lcv.md): genetic causality proportion requires sufficient heritability and LD-score inputs.
- [LHC-MR](../references/lhc-mr.md): bidirectional effects with heritability and confounding need study-specific formatting and reference preparation.
- Bidirectional MR: independently select and harmonize instruments in each direction, apply Steiger filtering and report both directions; reverse direction requires its own exposure instruments.
- Contamination mixture and MRMix: source describes these as robust alternatives; adapt the study-specific inputs and assess instrument assumptions. The reporting reference gives an interpretation framework.
- Assess NOME with I-squared; the source notes thresholds of 0.6 and 0.9 and a manual fallback if Isq is unavailable. See [failure modes](../references/failure-modes.md).

These are population-level research outputs, not individual treatment or diagnostic recommendations.

Related shelf skills: bio-causal-genomics-mendelian-randomization, bio-causal-genomics-genetic-correlation and bio-causal-genomics-colocalization-analysis.
