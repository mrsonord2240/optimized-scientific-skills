## Per-Format Failure Modes

### SMILES -- ambiguous aromaticity

**Trigger:** Input from non-RDKit source (ChemAxon, OpenEye, Daylight) round-tripping into RDKit.

**Mechanism:** RDKit perceives aromaticity on input. Aromatic flags from origin toolkit are overwritten.

**Symptom:** Fingerprints differ between toolkits for "identical" molecules; database joins by canonical SMILES miss records.

**Fix:** Always re-canonicalize within the analysis toolkit. For cross-toolkit identity, use InChIKey not canonical SMILES.

### SDF V2000 -- atom count >999

**Trigger:** Large molecules (peptides, oligonucleotides, dendrimers).

**Mechanism:** V2000 header uses fixed 3-character atom count field.

**Symptom:** Current RDKit writers automatically emit V3000 for oversized molecules, so the usual symptom is a downstream/legacy reader rejecting the resulting V3000 record. A workflow that explicitly writes V2000 cannot represent more than 999 atoms.

**Fix:** Treat V3000 as the required interchange format, and confirm every downstream reader accepts it. RDKit auto-detects V3000 on read; request it explicitly when the format choice must be visible and reproducible:

```python
writer = Chem.SDWriter('out.sdf')
writer.SetForceV3000(True)
writer.write(mol)
writer.close()
```

### SDF -- wedge bond orientation lost

**Trigger:** SDF written by tools that use parity flags only (older ISIS-Draw, some pipeline tools).

**Mechanism:** Parity alone is ambiguous without geometric coordinates; RDKit reads parity but cannot re-render wedges.

**Symptom:** Drawn molecule shows undefined stereo despite SDF carrying parity bits.

**Fix:** After read, `Chem.AssignStereochemistryFrom3D(mol)` if 3D coords present; otherwise stereo must be re-derived from SMILES with wedges.

### PDB ligand -- no bond orders

**Trigger:** Parsing ligand from PDB entry (e.g., extracting co-crystal ligand).

**Mechanism:** PDB stores only atoms + CONECT; bond orders inferred by RDKit's `AssignBondOrdersFromTemplate` which requires a template molecule.

**Symptom:** All bonds single; aromatic rings non-aromatic; valences wrong.

**Fix:** Use `AllChem.AssignBondOrdersFromTemplate(template, ligand)` where `template` is a SMILES-derived mol of the expected ligand structure. Or use the PDB Ligand Expo SDF.

### MOL2 -- SYBYL atom type dialect

**Trigger:** MOL2 produced by Corina, MOE, or Schrodinger.

**Mechanism:** SYBYL atom types are not perfectly standardized across vendors; RDKit's parser handles canonical SYBYL.

**Symptom:** Mol returns as `None` or with wrong atom types (`Cl` vs `Cl.O` peroxide-style).

**Fix:** Convert via Open Babel as intermediate: `obabel input.mol2 -O temp.sdf` then read SDF.

### Open Babel pybel -- import path

**Trigger:** Code written for Open Babel 2.x.

**Mechanism:** OB 3.x reorganized: `import pybel` no longer works.

**Symptom:** `ModuleNotFoundError: No module named 'pybel'`.

**Fix:** `from openbabel import pybel`.

## Charge Models on I/O

| Source | Charges in file | Use for |
|--------|-----------------|---------|
| Parsed SMILES | Atom-local formal charges; no partial-charge model | Storage, similarity, ML training |
| Parsed PDB | Atomic charges typically absent | Always re-assign for downstream |
| `obabel --partialcharge gasteiger` | Gasteiger-Marsili partial charges (empirical) | Workflows that explicitly require Gasteiger charges; Vina/Vinardo scoring itself does not require assigned atom charges |
| AM1-BCC (AmberTools antechamber) | Semi-empirical | MD, FEP setup |
| RESP (psi4, Gaussian) | Restrained fit to a quantum-mechanical ESP; protocol-specific | Force-field workflows parameterized for that RESP protocol |

The charge model **must** match the downstream method. Mixing AM1-BCC ligand charges with TIP3P water + AMBER protein is valid; Gasteiger charges are unsuitable for MD.

## Drawing for QC

Always draw a random subset of parsed molecules. Wrong stereo, missing rings, and broken aromaticity show immediately.

```python
from rdkit.Chem.Draw import rdMolDraw2D

def draw_grid(mols, fname, mols_per_row=5, sub_img_size=(250, 200)):
    from rdkit.Chem.Draw import MolsToGridImage
    img = MolsToGridImage(mols[:25], molsPerRow=mols_per_row, subImgSize=sub_img_size,
                          legends=[m.GetProp('_Name') if m.HasProp('_Name') else ''
                                   for m in mols[:25]])
    img.save(fname)
```

`MolsToGridImage` returns PIL image; for headless servers use `MolDraw2DCairo` directly.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `Chem.MolFromSmiles` returns None even with `sanitize=False` | Grammar error: invalid token, bad parentheses, or ring not closed | Correct or reject the source string; there is no molecule to inspect |
| `Chem.MolFromSmiles` returns None, but `sanitize=False` returns a molecule | Chemistry error: valence, kekulization, or aromaticity | Run `Chem.SanitizeMol(mol, catchErrors=True)` to identify the failing step |
| Empty/whitespace SMILES becomes a zero-atom `Mol` | Empty field is accepted as an empty graph | Reject blank text before parsing and reject `mol.GetNumAtoms() == 0` |
| Round-trip SMILES changes | Aromaticity perception drift | Always canonicalize within analysis toolkit |
| All bonds single in PDB ligand | PDB has no bond orders | `AllChem.AssignBondOrdersFromTemplate(template, mol)` |
| Stereo lost on SDF write | Stereo was absent, removed, or not represented by coordinates/bond directions | Verify assigned chiral tags and bond stereo before writing; preserve suitable 2D/3D coordinates and inspect the round trip |
| MOL2 parse returns None | RDKit MOL2 parser incomplete for vendor dialects | Convert via Open Babel intermediate |
| InChI differs for "same" molecule | Different tautomers, charges, or stereo | Use `/FixedH` to retain tautomer; compare without standardization |
| Fingerprints differ across toolkits | Aromaticity model difference | Use InChIKey for identity; re-canonicalize for similarity |
