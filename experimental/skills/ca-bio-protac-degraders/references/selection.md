## E3 Ligase Choice

| Recruited UPS component | Ligand series | Published design context | Limitations |
|-----------|---------------|---------|-------------|
| VHL | VL-269 (Gechijian et al. 2018) | Published VHL-recruiting degraders | Expression and productive geometry are system-dependent |
| CRBN (cereblon) | thalidomide, pomalidomide | Extensively used recruiter series | Neosubstrate liabilities depend on recruiter and context |
| IAP (XIAP, cIAP1) | SMAC-mimetic-derived recruiters | Published IAP-recruiting degraders | Target scope and cellular effects require validation |
| MDM2 | nutlin-derived recruiters | Published MDM2-recruiting degraders | Target diversity and pathway effects require validation |
| KEAP1 | KEAP1-directed recruiters | CUL3-KEAP1 recruitment studies | Specialized use and limited comparative validation |
| DCAF15 | Aryl sulfonamides such as E7820 | DDB1-CUL4 / DCAF15 systems | Molecular-glue and degrader mechanisms require careful distinction |
| RNF114 | Nimbolide, EN219 | Covalent RNF114 recruitment | Limited tooling |
| RNF4 | CCW16 | Covalent RNF4 recruitment | Limited tooling |
| UBE2D (E2, not E3) | EN450 | Covalent molecular-glue mechanism involving NFKB1 | Do not classify as an E3-ligase recruiter |

**Decision:** Select an E3 recruiter using evidence for ligand availability, target/E3 geometry, cellular expression, neosubstrate liabilities, and the intended biological system. CRBN and VHL are common starting points with extensive published examples, but neither is a universal first choice.

## Linker Design Principles

Linkers tune ternary complex geometry and stability. The ranges below are exploratory starting points, not validated acceptance criteria:

| Property | Range | Effect |
|----------|-------|--------|
| Linker length | Project-defined enumerated series | Critical; geometry-dependent |
| Linker rigidity | Flexible (PEG) vs rigid (piperazine, pyridine) | Changes the accessible conformational ensemble |
| Linker chemistry | PEG, alkyl, piperazine, triazole, ether, amide | PEG common; rigid for tighter binding |
| Click chemistry compatibility | Triazole-forming routes are one option | Requires route- and attachment-specific synthesis review |
| Molecular size and polarity | Measure across the designed series | Permeability, solubility, and exposure depend on the complete molecule and its conformations |

**Critical:** The "Goldilocks linker length" is target-specific. Too short can create a ternary clash; too long can impose an unfavorable entropic cost or permit unproductive geometries. Enumerate a series around the geometry supported by the binary structures rather than assuming a universal optimal range.

## Decision Tree by Scenario

| Goal | E3 / linker | Tools |
|------|-------------|-------|
| Initial degrader series | Compare supported recruiter and linker variants | PRosettaC for ternary hypotheses |
| Reduce recruiter-specific liabilities | Compare alternative E3 recruiters and linker geometries | Structural hypotheses + cellular selectivity validation |
| Target with prior recruiter-specific evidence | Reproduce the supported recruiter context, then vary deliberately | Match the published and intended biological systems |
| Targeted protein degradation program | Select E3 using geometry, expression, and liabilities | Structural and experimental validation track |
| Novel target without an established ternary model | Multiple E3 / linker variants | Combinatorial design + PRosettaC |
| Molecular glue (non-PROTAC) | Use a glue-specific discovery strategy | Distinct mechanism; do not treat as linker design |
| Characterize cooperativity | Structural hypotheses plus experiment | ITC or SPR/BLI with matched binary and ternary measurements |
| Cell-active candidate | Standard development | PK + degradation cellular assays |
