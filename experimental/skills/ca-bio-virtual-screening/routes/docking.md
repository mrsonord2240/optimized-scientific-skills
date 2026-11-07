# Vina docking

```bash
python scripts/dock_single.py receptor.pdbqt ligand.pdbqt --center 0 0 0 --size 20 20 20 --exhaustiveness 8 --n-poses 10 --seed 42 --out poses.pdbqt
```

Replace center and size with a validated box. The script prints only negative affinities but writes Vina's complete pose file; inspect and filter poses before downstream use. For a library, `scripts/virtual_screen.py` offers virtual_screen with a name-to-SMILES mapping and explicit center and box; see `references/vina.md`, `references/selection.md` and `references/failure-modes.md`. Benchmark self- and cross-docking controls.

Done when poses, scores, settings and failures are retained for review.
