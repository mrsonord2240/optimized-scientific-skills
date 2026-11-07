# ETKDGv3 ensemble

```python
from rdkit import Chem
from scripts.gen_conformers import gen_conformer_ensemble, prune_rmsd, filter_energy_window
smiles = open('input.smi').readline().strip()
mol, data = gen_conformer_ensemble(smiles, n_conf=20, seed=42)
kept = filter_energy_window(prune_rmsd(mol, data, rmsd_cutoff=0.5), window_kcal=10.0)
writer = Chem.SDWriter('conformers.sdf')
for cid, energy in kept:
    mol.SetProp('energy_kcal_mol', str(energy))
    writer.write(mol, confId=cid)
writer.close()
```

The helper is `scripts/gen_conformers.py`; it checks MMFF94s parameters, falls back to parameter-checked UFF, and refuses nonconvergence. These pruning parameters are starting heuristics; validate coverage and the energy model for the task. For method selection see `references/selection.md`; see `references/etkdg.md` and `references/troubleshooting.md`.

Done when conformers and settings are saved and sampling adequacy assessed.
