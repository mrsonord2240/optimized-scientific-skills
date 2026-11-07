# Ternary prediction tools: how to invoke each, and PRosettaC vs AlphaFold3

## How to invoke each tool

| Tool | How to invoke | Turnaround | What you get back |
|------|----------------|------------|--------------------|
| PRosettaC | Web submission at prosettac.weizmann.ac.il (free registration); upload target PDB + target-ligand pose, E3 PDB + E3-ligand pose, and the PROTAC SMILES | Hours to ~1 day (queued Rosetta job) | Ranked/clustered ternary poses + interface scores; download the top pose's exit-vector distance to sanity-check candidate linkers against `ternary_geometry_screen.py` |
| AlphaFold3 | AlphaFold Server (web, request access) for a quick check, or a licensed local install (request model weights from Google DeepMind) for batch/automated use; submit full protein sequences + PROTAC SMILES as one unrestrained complex | Minutes (server) to GPU-hours (local) | A single predicted complex + per-residue confidence (pLDDT/PAE); no arbitrary distance restraint, so compare the predicted interface against known ternary structures rather than trusting confidence alone |
| Boltz-1 / Boltz-2 | Open-source, pip-installable (`pip install boltz`), but needs a local CUDA GPU and a multi-GB weights download on first run; `boltz predict input.yaml` per the project's README | GPU-minutes locally | A predicted complex structure, same caveats as AlphaFold3 (limited PROTAC-specific validation) |
| DeepTernary | Clone the research repo from GitHub and follow its own install/checkpoint instructions; no hosted API, no pip package | GPU-minutes locally, plus one-time setup | SE(3)-equivariant ternary pose; validate against known structures before using it to rank |
| HADDOCK | HADDOCK3 local install (compiles against CNS, which is itself registration-gated academically), or the HADDOCK web portal at bonvinlab.org (free academic account, upload structures + manual restraints) | Local: setup-heavy. Web: minutes to hours queued | Restraint-guided docking poses; you must supply the restraints yourself |

## Reconciliation: PRosettaC vs AlphaFold3

| Aspect | PRosettaC | AlphaFold3 |
|--------|-----------|------------|
| Approach | PROTAC-specific Rosetta sampling | Unrestrained foundation-model prediction |
| Accuracy | Higher average DockQ in one 36-structure comparison, but only 25 complexes were modeled and most predictions were low quality | Limited PROTAC-specific validation |
| Speed | Measure for the installed workflow and hardware | Measure for the selected service or local hardware |
| Access | Web service | AlphaFold Server or local installation, subject to their terms and limits |
| Restraints | Method-specific setup | No arbitrary user distance restraints |
| Decision | Use as a PROTAC-specific structural hypothesis | Use as an independently benchmarked structural hypothesis |

Use PRosettaC or another benchmarked structural method to generate hypotheses, then measure ternary binding/cooperativity and cellular degradation experimentally. Do not treat any one modeling method as a validated universal ranker.
