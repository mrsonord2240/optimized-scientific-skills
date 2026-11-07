## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Vina segfault | PDBQT corrupted (atom names) | Re-prep with meeko |
| GNINA hangs | GPU OOM | Reduce concurrent work and, if fewer output poses are acceptable, use `--num_modes 5` |
| All affinities very poor (-3 to -5) | Wrong protonation; ligand too large for box | Re-check pKa; expand box |
| Identical affinity across ligands | Receptor grid not computed | Call `v.compute_vina_maps()` before dock |
| PoseBusters passes only ~3/12 on a docked charged ligand; `Explicit valence ... is greater than permitted` | `obabel` PDBQT -> SDF dropped bond orders/charges | Rebuild with meeko `RDKitMolCreate.from_pdbqt_mol` (see Handoff caveat) |
| `prank predict` exits with no `pockets/` directory and no error | `prank.bat` returns immediately when run through a shell wrapper on Windows | Call the jar directly: `java -cp ".../bin/p2rank.jar;.../bin/lib/*" cz.siret.prank.program.Main predict ...` (see Binding Site Detection) |
| Pose poses make no sense | Receptor and ligand in different frames | Ensure same coordinate origin |
| Metal-coordination pose is wrong | The selected scoring/preparation protocol lacks a validated model for that metal geometry | Use a metal-specific validated workflow; the Vina executable can use AutoDock4Zn maps with `--scoring ad4` for zinc, while other metals require separately supported parameters/protocols |
| GPU mode slow | Vina is CPU-only; only GNINA is GPU | Use GNINA for GPU; if using a third-party GPU port of Vina, benchmark it on the same hardware, target, library tranche, and search settings before adopting it |
| `mk_prepare_receptor.py: command not found` | meeko's pip-installed console-script has no `.py` suffix | Call `mk_prepare_receptor` (no `.py`), as `scripts/prepare_receptor.py` does |
| `ValueError: invalid literal for int() with base 10: '184A'` from `mk_prepare_receptor --read_pqr` | Any residue with a PDB insertion code (e.g. chymotrypsin-numbered serine proteases such as trypsin, elastase): meeko's PQR reader assumes an all-integer residue-number column | Use `pdb2pqr --pdb-output` + `mk_prepare_receptor --read_pdb` instead of `--read_pqr` |
| `pip install vina` fails with "Boost library location was not found" | No Windows wheel for the `vina` PyPI package | Use the Vina CLI via `subprocess` instead of `from vina import Vina` |
