# Receptor and ligand preparation

```bash
python scripts/prepare_receptor.py repaired.pdb receptor.pdbqt --ph 7.4
```

Review retained waters, cofactors, metals and missing atoms first. For ligand preparation, choose a documented protomer then use `scripts/virtual_screen.py` as shown in `references/ligand.md`. Receptor preparation writes protonated intermediate files beside the input; see `references/receptor.md`.

Done when prepared structures and protonation assumptions are recorded.
