## Example: KRAS G12C Inhibitor Design Workflow

**Goal:** Decorate a co-crystal scaffold with a cysteine-targeting warhead and rank candidates by covalent efficiency.

**Approach:** Load scaffold SMILES, enumerate acrylamide-bearing analogs, filter by reactivity selectivity, dock under covalent constraint, and rank by kinact/Ki surrogates.

```python
from rdkit import Chem

# Step 1: scaffold from co-crystal (4LRW or AMG510)
scaffold_smi = 'c1ccc(C(=O)NCC)cc1'  # generic valid scaffold for code illustration
scaffold = Chem.MolFromSmiles(scaffold_smi)

# Step 2: enumerate analogs with acrylamide warhead
def add_acrylamide(scaffold, attachment_atom_idx):
    """Project hook: attach a mapped acrylamide with an audited reaction."""
    raise NotImplementedError(
        'Provide a project-specific mapped reaction and validate atom mapping, '
        'valence, regioisomer identity, and product sanitization.'
    )

# Step 3: filter for reactive group selectivity
# Step 4: dock with DOCKovalent / GOLD covalent / HCovDock
# Step 5: rank by kinact/Ki surrogate (compute reactive Michael acceptor reactivity)
```

## Reactivity Surrogates (computed without experiment)

For ranking warheads without wet-lab data:

| Descriptor | Use case |
|------------|----------|
| LUMO energy (DFT) | Michael acceptor reactivity (lower LUMO = more reactive) |
| Electrophile partial charge | SN2 reactivity |
| RDKit `rdMolDescriptors.CalcLabuteASA` | Steric accessibility |
| Experimentally supported binding pose or validated docking model | Geometric fit to reactive residue |

**Goal:** Record alpha-carbon substitution as a structural feature for an acrylamide series.

**Approach:** Parse the SMILES, locate the acrylamide substructure, and count neighbors on the alpha carbon outside the matched warhead. This count is not a LUMO estimate or a stand-alone reactivity prediction; substituent electronics and the rest of the molecule must be considered, and reactivity should be measured.

Run `python scripts/acrylamide_alpha_substitution.py '<smiles>' ...` (or import `acrylamide_alpha_substitution_count` from it); it prints the count, or `None` for an unparsable SMILES or no acrylamide.

For a reactivity model, use experimentally measured rates or a validated quantum-chemical workflow; a single frontier-orbital energy is not sufficient on its own.
