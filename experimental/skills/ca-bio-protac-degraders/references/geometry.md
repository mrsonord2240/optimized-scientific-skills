## Ternary Complex Modeling Workflow

**Goal:** Predict 3D structure of target-PROTAC-E3 ternary complex.

**What this Skill runs locally (linker pre-filter only, not a structure prediction):**
1. Start with binary co-crystals: target + target-ligand pose; E3 + E3-ligand pose
2. Measure the target-ligand / E3-ligand exit-vector distance yourself with
   `attachment_distance()` (below) once both are placed in a shared coordinate frame
3. Enumerate candidate linkers (`scripts/protac_enumerate.py`) and screen them against that
   required span with `scripts/ternary_geometry_screen.py`'s `screen_linker_library()` -- it
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

Runnable: `scripts/ternary_geometry_screen.py`'s `linker_reach()` embeds standalone-linker
conformers (RDKit ETKDG + MMFF) and reports the min/mean/max attachment-point distance each
named linker can achieve; `screen_linker_library()` compares that range against a required span
(from `attachment_distance()` on your own structures) with a tolerance. This has no knowledge of
the target/E3 protein context -- it bounds linker reach only, not ternary-complex feasibility.

## Generative Linker Design

REINVENT 4 can generate linkers, but it does not provide the `ternary_score` / `deepternary` interface shown in some informal examples. Export generated candidates, run an installed and validated ternary-prediction workflow separately, and then join the structural scores back to the candidates. Do not assume DeepTernary is a web API or a built-in REINVENT scoring component.
