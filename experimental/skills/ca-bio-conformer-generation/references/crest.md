## CREST + GFN2-xTB for High-Quality Sampling

CREST (Pracht et al. 2024) performs iterative meta-dynamics + GFN2-xTB optimization for conformer sampling.

**Platform:** CREST is Linux/macOS only -- no Windows build (see Installation). On native Windows, run it under WSL, or use the RDKit macrocycle-aware fallback above plus standalone `xtb --opt`.

**Goal:** Sample high-quality conformer ensembles for macrocycles, peptides, or molecules where ETKDGv3 + MMFF94 is inadequate.

**Approach:** Start from an RDKit-generated MMFF94-relaxed conformer, write to XYZ, and run CREST with GFN2-xTB driver to perform iterative meta-dynamics + reoptimization.

```bash
xtb mol.xyz --opt extreme
crest xtbopt.xyz --gfn2 --T 12 --ewin 6
```

**`--gfn2`**: use GFN2-xTB; validate its coverage and energy ordering for the chemistry of interest.
**`--gfn-ff`**: use GFN-FF; it can reduce computational cost, but benchmark its coverage and energy ordering for the intended molecules.
**`-ewin 6`**: 6 kcal/mol energy window above global min.
**`-T 12`**: use 12 CPU threads.

Output: `crest_conformers.xyz` with sampled ensemble.

**Workflow:** Start from RDKit ETKDGv3 + MMFF94 (cheap initial structure) -> save as XYZ -> CREST refinement.

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from pathlib import Path

def crest_workflow(smiles, out_dir='crest_out'):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f'Invalid SMILES: {smiles!r}')
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    if AllChem.EmbedMolecule(mol, params) == -1:
        raise RuntimeError(f'Initial 3D embedding failed for {smiles!r}')
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError(f'MMFF parameters are unavailable for {smiles!r}')
    status = AllChem.MMFFOptimizeMolecule(
        mol, mmffVariant='MMFF94s', maxIters=1000
    )
    if status != 0:
        raise RuntimeError(f'Initial MMFF optimization did not converge (status {status})')

    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    input_path = out_dir / 'input.xyz'
    input_path.write_text(Chem.MolToXYZBlock(mol))
    # Because cwd is out_dir, pass the coordinate filename rather than out_dir/input.xyz.
    subprocess.run(['crest', input_path.name, '--gfn2', '-T', '12'],
                   cwd=out_dir, check=True)
    output_path = out_dir / 'crest_conformers.xyz'
    if not output_path.exists():
        raise FileNotFoundError(f'CREST did not produce {output_path}')
    return output_path
```
