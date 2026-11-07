# Relative ligand strain

```python
from scripts.validate_poses import ligand_strain_mmff
ligand_strain_mmff('poses.sdf', n_ref_conf=20).to_csv('strain.csv', index=False)
```

This compares MMFF energy with a sampled reference minimum after relaxing only added hydrogens at fixed heavy atoms. It is not solution free energy. Review missing parameters and sampling coverage; use only a project-validated cutoff. See `references/strain.md`.

Done when strain and unscored poses are reported.
