## Vina Docking (Single Ligand)

```bash
python scripts/dock_single.py receptor.pdbqt ligand.pdbqt --center X Y Z --size X Y Z \
       --exhaustiveness 8 --n-poses 10 --seed 42 --out poses.pdbqt
```

`scripts/dock_single.py` uses the Vina Python API (Vina 1.2+) when `from vina import Vina` imports, and otherwise runs the Vina CLI (`--vina-exe` if `vina` is not on `PATH`; Windows has no `vina` wheel, see Version Compatibility). Both paths pass the same seed. It writes the poses to `--out` and prints one affinity per mode with `affinity >= 0` modes dropped (Vina occasionally emits a physically nonsensical positive-energy mode, e.g. +68 kcal/mol; see "Sanity-filter reported poses").

**Exhaustiveness:** `8` is the Vina default. Increasing it increases search effort, but runtime and pose recovery depend on hardware, ligand flexibility, box size, and software version. Benchmark settings such as 8, 16, 32, and 64 on target-relevant controls instead of assigning universal timing or quality labels.

**Seed and reproducibility:** `scripts/dock_single.py` and `scripts/virtual_screen.py` accept a `seed` (CLI: `--seed`). Top-1 affinity is empirically stable run-to-run without a fixed seed, but poses ranked 2+ reorder between runs on the same input. Set and record a seed (default `42` above) whenever the top-N poses -- not only the single best -- will be reported or compared.

**Sanity-filter reported poses:** Vina's raw mode list can include a physically nonsensical outlier (a positive-energy mode was observed among 9 returned modes in testing). Filter or flag `affinity >= 0` poses before reporting or ranking; do not assume every mode Vina returns is a plausible binder.

Vina's `rmsd_lb` and `rmsd_ub` are lower and upper heavy-atom RMSD bounds between a reported mode and the best-scoring mode; the bounds differ in how symmetry-equivalent atoms are handled. They are not pose-versus-experimental-reference RMSDs. Use an external symmetry-aware RMSD to a reference pose for accuracy QC.
