## Receptor Preparation

**Goal:** Convert a protein PDB into a docking-ready format with correct protonation, missing atoms, and removed waters.

**Approach:** Decide which ligands, cofactors, metals, and structural waters to retain -> fill missing heavy atoms with a structure-repair tool such as PDBFixer -> use PROPKA/PDB2PQR plus manual review to assign pH-dependent protonation -> assign the charge model required by the docking workflow -> prepare receptor PDBQT with a documented AutoDock-compatible tool.

```bash
python scripts/prepare_receptor.py repaired.pdb receptor.pdbqt --ph 7.4
```

`scripts/prepare_receptor.py` runs `pdb2pqr --ff=AMBER --with-ph=<pH> --pdb-output` then `mk_prepare_receptor --read_pdb ... -p`. Decide which waters, cofactors and metals to retain before calling it. It needs `pdb2pqr` and `mk_prepare_receptor` on `PATH`; it writes the protonated `<base>_pH<pH>.pdb`/`.pqr` beside the input. It hands meeko the protonated PDB rather than the `.pqr` (see the insertion-code pitfall below).

**Common pitfall:** Forgetting to add hydrogens at protein pH (7.4) but using pH 7.0 ligand charges. Hist mistakenly protonated. Use PROPKA + manual review of catalytic residues.

**Common pitfall:** Feeding pdb2pqr's default `.pqr` output straight into `mk_prepare_receptor --read_pqr`. meeko 0.8.0's PQR reader assumes an all-integer residue-number column and raises `ValueError: invalid literal for int() with base 10` on any residue with a PDB insertion code (e.g. `184A`). Any PDB deposition with insertion-code residues triggers it, not one protein family: chymotrypsin-numbered serine proteases (trypsin 3PTB: 184A/188A/221A; elastase 1EAI) are common examples and hit it on real structures, not just edge cases. Use `pdb2pqr --pdb-output` and `mk_prepare_receptor --read_pdb` as above; verified on PDB 3PTB (trypsin, insertion-code residues 184A/188A/221A).
