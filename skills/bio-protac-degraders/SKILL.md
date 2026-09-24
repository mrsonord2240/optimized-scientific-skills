---
name: bio-protac-degraders
category: Data Analysis
description: Designs PROTACs, molecular glues, and bivalent degraders with explicit handling of E3 ligase choice (VHL, CRBN, IAP, MDM2, KEAP1), linker design (length, composition, rigidity), ternary complex prediction (PRosettaC, DeepTernary, AlphaFold3), cooperativity (alpha), DC50 / Dmax characterization, hook effect, and prediction-experiment reconciliation. Use when designing targeted protein degraders, planning linker SAR, predicting ternary complex stability, or building generative degrader workflows.
tool_type: python
primary_tool: PRosettaC
license: MIT
author: GPTomics
---

## Version Compatibility

The three shipped scripts (`examples/protac_enumerate.py`, `examples/ternary_geometry_screen.py`,
`examples/cooperativity_dc50.py`) need only RDKit, numpy, and scipy -- checked on RDKit 2026.03.6,
numpy 2.x, scipy 1.18. RDKit 2024.09+ should also work via the introspect-and-adapt pattern below.

PRosettaC, DeepTernary, AlphaFold3, Boltz-1/2, and HADDOCK are external services or
GPU/licence-gated tools this Skill does not install or run locally -- see "Ternary Complex
Prediction Tools" below for how to invoke each one directly, and OpenMM 8.1+ if you refine a
ternary pose returned by one of them.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# PROTAC and Bivalent Degrader Design

Design bifunctional molecules (PROTACs) that recruit an E3 ubiquitin ligase to a target protein, inducing target ubiquitination and proteasomal degradation. PROTACs differ from traditional drugs: a productive **ternary complex** (target + PROTAC + E3) is required, not just target binding. The modality has produced clinical programs, but their development and regulatory status changes rapidly and must be checked from current sources. PROTAC design balances **target ligand binding**, **E3 ligand binding**, **linker geometry** (length, rigidity, chemistry), **cooperativity**, **dose-dependent ternary-complex formation**, and **cell permeability**. Negative cooperativity and the high-concentration hook effect are distinct phenomena, although cooperativity can influence the dose-response profile.

For target ligand design, see `chemoinformatics/virtual-screening` and `chemoinformatics/admet-prediction`. For linker-only enumeration, see `chemoinformatics/reaction-enumeration`. For generative linker design, see `chemoinformatics/generative-design`.

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

## Ternary Complex Prediction Tools

| Tool | Approach | Strength | Fails when |
|------|----------|----------|------------|
| PRosettaC | Constrained PatchDock, RosettaDock refinement, PROTAC conformer generation, repacking, and clustering | PROTAC-specific published workflow | Performance varies by complex; Rosetta/service requirements |
| DeepTernary | Equivariant deep learning | Fast; SE(3) | OOD chemistry |
| AlphaFold3 | Unrestrained whole-complex prediction | Accepts proteins and ligands | No arbitrary user distance-restraint interface; benchmark PROTAC use |
| Boltz-1 / Boltz-2 | Unrestrained whole-complex prediction | Open local models | Limited PROTAC-specific validation |
| HADDOCK | Information-driven, restraint-guided docking | Mature integrative docking framework | Manual restraint specification |

**Decision:** Use a PROTAC-specific method such as **PRosettaC** for first-pass ternary modeling. AlphaFold3 or Boltz can provide unrestrained whole-complex predictions, but should be benchmarked on relevant ternary complexes. **DeepTernary** is released research code rather than a hosted API; validate it against relevant structures before prospective ranking.

**None of these tools run locally through this Skill.** All are external submissions or gated
local installs -- this Skill's only local, runnable content for ternary hypotheses is
`examples/ternary_geometry_screen.py`, a linker-reach pre-filter (see "Ternary Complex Modeling
Workflow" below); it does not predict a structure or an interface score. To actually run one of
the methods above, read `references/ternary-prediction-tools.md` (per-tool URL, account, turnaround, output;
also a PRosettaC vs AlphaFold3 comparison).

Do not fabricate a numeric score from any of these tools when they have not actually been run --
state plainly that the step requires the external submission above.

## Cooperativity (Alpha)

Cooperativity quantifies how the ternary complex stabilizes (or destabilizes) the binary binding:

```
alpha = (Kd_binary,target) / (Kd_ternary,target)
```

- alpha > 1: positive cooperativity (ternary stronger than binary)
- alpha = 1: no cooperativity (independent binding)
- alpha < 1: negative cooperativity (mutual destabilization)

Positive cooperativity can favor ternary-complex formation, but the preferred alpha is system- and assay-dependent and does not alone establish degradation efficacy. Alpha must be measured from binary and ternary binding experiments; a predicted structure does not directly provide it.

Measure with ITC (isothermal titration calorimetry) or SPR/BLI titrations of binary vs ternary.

Runnable: `examples/cooperativity_dc50.py`'s `cooperativity_alpha(kd_binary, kd_ternary)` applies
this formula to your measured Kd pair and labels positive/negative/no cooperativity.

## DC50 / Dmax Characterization

In cellular assays:
- **DC50**: PROTAC concentration for 50% degradation (analogous to IC50)
- **Dmax**: maximum fraction degraded at any concentration

| Property | What to report | Interpretation |
|----------|----------------|----------------|
| DC50 | Concentration producing 50% of the assay's fitted maximal degradation | Compare only across matched assay conditions; no universal clinical cutoff |
| Dmax | Maximum observed or fitted degradation and uncertainty | Required depletion is target- and phenotype-dependent |
| Hook effect | Full concentration-response range and concentration of any downturn | A high-concentration effect; its location is system- and assay-dependent |
| Cooperativity | Alpha from matched binary and ternary binding experiments | Distinct from the hook effect and insufficient by itself to predict degradation |

**Hook effect**: at high PROTAC concentrations, binary complexes (PROTAC-target alone, PROTAC-E3 alone) dominate, and ternary complex formation drops. Dose-response curves are bell-shaped.

Runnable: `examples/cooperativity_dc50.py` fits DC50/Dmax with a 3-parameter Hill curve, flags a
hook effect from a >15-percentage-point downturn after the peak, and -- when a hook is flagged --
restricts the DC50/Dmax fit to the ascending arm rather than fitting a monotonic sigmoid to
non-monotonic data. When the peak sits under 10x the fitted DC50 (hook onset close to DC50), the
ascending arm never reaches its plateau and Dmax is underestimated (15 pp on a noise-free
synthetic curve); `fit_dc50()` then returns `plateau_reached=False` and a `caveat` -- report Dmax as
a lower bound. It ships a seeded synthetic dose-response curve and checks the fit recovers
the curve's own planted DC50/Dmax within tolerance, and that a narrow-separation hook is caveated.

## Ternary Complex Modeling Workflow

**Goal:** Predict 3D structure of target-PROTAC-E3 ternary complex.

**What this Skill runs locally (linker pre-filter only, not a structure prediction):**
1. Start with binary co-crystals: target + target-ligand pose; E3 + E3-ligand pose
2. Measure the target-ligand / E3-ligand exit-vector distance yourself with
   `attachment_distance()` (below) once both are placed in a shared coordinate frame
3. Enumerate candidate linkers (`examples/protac_enumerate.py`) and screen them against that
   required span with `examples/ternary_geometry_screen.py`'s `screen_linker_library()` -- it
   embeds 3D conformers of each linker and reports whether its sampled reach range can plausibly
   span the requirement, as a cheap filter before spending an external ternary submission
4. **Structure prediction/refinement itself is external** -- submit the surviving candidates to
   PRosettaC, AlphaFold3, Boltz, DeepTernary, or HADDOCK per `references/ternary-prediction-tools.md`; none of
   that happens locally

For a production workflow, use PRosettaC or provide both proteins and the complete PROTAC as components of an unrestrained AlphaFold3 input. AlphaFold3 does not expose arbitrary chain-chain distance restraints; compare predicted interfaces and confidence with known complexes or a PROTAC-specific method.

## Linker Geometry Assessment

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def attachment_distance(target_ligand, e3_ligand,
                        target_attachment_idx, e3_attachment_idx):
    """
    Measure an attachment-point distance after both ligands have been placed in
    the same ternary-complex coordinate frame.
    """
    p1 = target_ligand.GetConformer().GetAtomPosition(target_attachment_idx)
    p2 = e3_ligand.GetConformer().GetAtomPosition(e3_attachment_idx)
    return p1.Distance(p2)
```

An attachment-point distance does not map uniquely to a linker atom count: bond geometry, rigidity, branching, solvation, and the relative protein orientation all matter. Enumerate chemically synthesizable linker candidates, sample their conformers in the ternary geometry, and retain candidates that can connect without severe strain or clashes.

Runnable: `examples/ternary_geometry_screen.py`'s `linker_reach()` embeds standalone-linker
conformers (RDKit ETKDG + MMFF) and reports the min/mean/max attachment-point distance each
named linker can achieve; `screen_linker_library()` compares that range against a required span
(from `attachment_distance()` on your own structures) with a tolerance. This has no knowledge of
the target/E3 protein context -- it bounds linker reach only, not ternary-complex feasibility.

## Generative Linker Design

REINVENT 4 can generate linkers, but it does not provide the `ternary_score` / `deepternary` interface shown in some informal examples. Export generated candidates, run an installed and validated ternary-prediction workflow separately, and then join the structural scores back to the candidates. Do not assume DeepTernary is a web API or a built-in REINVENT scoring component.

## Per-Tool Failure Modes

### PRosettaC -- inaccessible E3 in selected ligase

**Trigger:** Target's known binding mode incompatible with E3 ligase orientation.

**Mechanism:** The selected binary poses, exit vectors, linker conformations, or protein orientation may not support a compatible ternary geometry.

**Symptom:** Low ternary complex scores; high RMSD across replicates.

**Fix:** Try a different E3 such as CRBN or VHL and compare against experimentally resolved ternary complexes with compatible exit-vector geometry.

### DeepTernary -- novel chemotype

**Trigger:** Target ligand or E3 ligand outside training distribution.

**Mechanism:** A ligand, linker, target, or E3 outside the method's validated domain may require extrapolation.

**Symptom:** Predicted ternary complex unrealistic.

**Fix:** Compare with an independently configured structural method and relevant known complexes; validate prospective ranking experimentally.

### Hook effect at high PROTAC concentration

**Trigger:** PROTAC concentration becomes high enough that separate target-PROTAC and E3-PROTAC binary complexes compete with productive ternary-complex formation.

**Mechanism:** Saturation by binary complexes reduces the population of productive ternary complex. Negative cooperativity can worsen ternary formation but is not the definition of the hook effect.

**Symptom:** Degradation increases and then decreases across a sufficiently broad concentration-response experiment.

**Fix:** Confirm the downturn experimentally over a broad dose range, then optimize ternary-complex geometry, cooperativity, exposure, and dosing without assuming that linker shortening alone will solve it.

### Insufficient cell permeability

**Trigger:** Measured permeability or cellular exposure is poor relative to biochemical activity.

**Mechanism:** Size, exposed polarity, conformation, ionization, or efflux may limit intracellular exposure.

**Symptom:** Cellular degradation potency is substantially worse than biochemical ternary-complex or binding measurements.

**Fix:** Optimize linker and exposed polarity using measured permeability, solubility, and intracellular exposure across the series. Do not impose a universal MW or TPSA cutoff.

### E3-target distance miscalculation

**Trigger:** Computing linker length from binary models without ternary refinement.

**Mechanism:** A distance from separately positioned binary structures does not determine the accessible ternary geometry or linker conformational ensemble.

**Symptom:** PROTACs synthesized at wrong linker length; no degradation.

**Fix:** Use a ternary structural hypothesis to define a chemically diverse linker series, then compare conformational feasibility and experimental degradation across that series.

### Molecular glue vs PROTAC confusion

**Trigger:** Designing as PROTAC when target lacks defined ligand.

**Mechanism:** A molecular glue stabilizes or induces a protein-protein interaction without the two-ligand-plus-linker architecture assumed by a PROTAC workflow.

**Symptom:** Design too rigid; no degradation despite ternary prediction.

**Fix:** For targets without known ligand, consider molecular glue discovery instead.


## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| PRosettaC fails to converge | Input geometry, sampling, or service/configuration problem | Inspect inputs and logs; compare justified recruiter/linker hypotheses |
| DeepTernary returns clashing pose | Prediction outside a validated domain or incorrect interface | Inspect confidence and clashes; compare an independent method and known structures |
| AlphaFold3 ternary unrealistic | Unrestrained prediction has a low-confidence or incorrect interface | Inspect confidence and compare with PRosettaC or known ternary structures |
| Cellular phenotype disagrees with target-degradation assays | Exposure, off-target degradation, assay timing, or pathway effects | Measure target engagement/degradation and use proteome-wide selectivity assays where appropriate |
| Degradation decreases at high PROTAC concentration | Hook effect from competing binary complexes | Confirm with a broad dose range; optimize ternary geometry and exposure |
| Synthesis is impractical | Proposed connectivity lacks a credible route | Obtain medicinal-chemistry review and redesign attachment chemistry or linker |
| Poor permeability or intracellular exposure | Size, exposed polarity, conformation, or efflux | Measure the bottleneck and optimize the series; avoid a universal size cutoff |
| `protac_enumerate.build_protac()` raises "Attachment dummy bonds must be single bonds" | An explicit `[*]` exit-vector dummy is attached via a double bond or other non-single bond | Use a single-bond attachment point; only relax this if the intended connection chemistry has its own, separately validated bond-order rule |

## Reference Files

| File | Read when |
|------|-----------|
| `references/ternary-prediction-tools.md` | Submitting a ternary-complex job (PRosettaC, AlphaFold3, Boltz, DeepTernary, HADDOCK), or choosing between PRosettaC and AlphaFold3 |

## References

- Békés M, Langley DR, Crews CM. "PROTAC targeted protein degraders: the past is prologue." *Nat. Rev. Drug Discov.* 21:181–200 (2022). DOI: 10.1038/s41573-021-00371-6.
- Drummond ML, Williams CI. "In Silico Modeling of PROTAC-Mediated Ternary Complexes: Validation and Application." *J. Chem. Inf. Model.* 59:1634–1644 (2019). DOI: 10.1021/acs.jcim.8b00992.
- Schapira M, Calabrese MF, Bullock AN, Crews CM. "Targeted protein degradation: expanding the toolbox." *Nat. Rev. Drug Discov.* 18:949–963 (2019). DOI: 10.1038/s41573-019-0047-y.
- Gechijian LN et al. "Functional TRIM24 degrader via conjugation of ineffectual bromodomain and VHL ligands." *Nat. Chem. Biol.* 14:405–412 (2018). DOI: 10.1038/s41589-018-0010-y.
- Zaidman D, Prilusky J, London N. "PRosettaC: Rosetta Based Modeling of PROTAC Mediated Ternary Complexes." *J. Chem. Inf. Model.* 60:4894–4903 (2020). DOI: 10.1021/acs.jcim.0c00589.
- Xue F, Zhang M, Li S et al. "SE(3)-equivariant ternary complex prediction towards target protein degradation." *Nat. Commun.* 16:5514 (2025). DOI: 10.1038/s41467-025-61272-5.
- Schulz JM, Schürer SI, Reynolds RC, Schürer SC. "PRosettaC outperforms AlphaFold3 for modeling PROTAC ternary complexes." *Sci. Rep.* 15:37620 (2025). DOI: 10.1038/s41598-025-21502-8.
- Bondeson DP et al. "Catalytic in vivo protein knockdown by small-molecule PROTACs." *Nat. Chem. Biol.* 11:611–617 (2015). DOI: 10.1038/nchembio.1858.
- Ward CC et al. "Covalent Ligand Screening Uncovers a RNF4 E3 Ligase Recruiter for Targeted Protein Degradation Applications." *ACS Chem. Biol.* 14:2430–2440 (2019). DOI: 10.1021/acschembio.8b01083.
- Luo M et al. "Chemoproteomics-enabled discovery of covalent RNF114-based degraders that mimic natural product function." *Cell Chem. Biol.* 28:559–566.e15 (2021). DOI: 10.1016/j.chembiol.2021.01.005.
- King EA et al. "Chemoproteomics-enabled discovery of a covalent molecular glue degrader targeting NF-kappaB." *Cell Chem. Biol.* 30:394–402.e9 (2023). DOI: 10.1016/j.chembiol.2023.02.008.
- Abramson J et al. "Accurate structure prediction of biomolecular interactions with AlphaFold 3." *Nature* 630:493–500 (2024). DOI: 10.1038/s41586-024-07487-w.
- REINVENT 4 official repository and installation documentation: https://github.com/MolecularAI/REINVENT4.
- HADDOCK3 official documentation: https://www.bonvinlab.org/haddock3/.
- Boltz official repository and documentation: https://github.com/jwohlwend/boltz.

## Related Skills

- chemoinformatics/molecular-io - Parse linker and ligand SMILES
- chemoinformatics/reaction-enumeration - Linker enumeration combinatorial
- chemoinformatics/generative-design - REINVENT linker mode
- chemoinformatics/conformer-generation - Ternary conformer sampling
- chemoinformatics/virtual-screening - Validate target ligand binding
- chemoinformatics/free-energy-calculations - Ternary ABFE / cooperativity
- chemoinformatics/admet-prediction - PROTAC ADMET specific challenges
- structural-biology/structure-io - PDB / mmCIF for ternary complex
