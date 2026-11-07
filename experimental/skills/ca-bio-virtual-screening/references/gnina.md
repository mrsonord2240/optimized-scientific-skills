# GNINA reference (moved from SKILL.md)

**Verified 2026-09-21** on GNINA 1.3.3 (`gnina.cuda12.8.static`), run through WSL2 with real GPU
passthrough (this machine has an RTX 5070 Ti) — no Docker needed, contrary to an earlier pass's
assumption. Docked benzamidine against the 3PTB fixture at this Skill's own recorded box (center
-1.52/14.47/17.47, size 20/20/20 A, `--cnn_scoring rescore`): base affinity -6.01 kcal/mol, matching
the Vina-only run's -5.978 on the same fixture, with a real CNN pose score (0.9646) and CNN affinity
(3.593) confirming the neural net actually scored on GPU. See the env's `TOOLS.md` for the install
(`cuda-toolkit=12.8` as one package, not individual libraries).

## GNINA with CNN Scoring (modern default)

```bash
gnina -r receptor.pdb -l ligand.sdf \
      --autobox_ligand reference_ligand.sdf \
      --cnn_scoring rescore \
      -o poses.sdf.gz \
      --num_modes 9 --exhaustiveness 8
```

`--cnn_scoring`:
- `none`: no CNN; use the selected empirical scoring function throughout
- `rescore` (default): use empirical scoring during the search, then CNN-rerank the final poses; least computationally expensive CNN option
- `refinement`: use the CNN to refine poses after Monte Carlo chains and to rank the final poses; approximately 10 times slower than `rescore` on a GPU in the official documentation
- `metrorescore`: use CNN scoring in the Metropolis search and rescore the resulting poses
- `metrorefine`: use CNN scoring in the Metropolis search and refine the resulting poses
- `all`: use the CNN scoring function throughout; the official documentation describes this as extremely computationally intensive and not recommended

The six choices above are from GNINA 1.3. Earlier releases expose a smaller set; check `gnina --help` for the installed executable rather than assuming every mode is available.

`--autobox_ligand`: define box from reference ligand SDF/PDB. Otherwise specify `--center_x/y/z` + `--size_x/y/z`.

**Critical:** GNINA distributions include multiple named CNN models/ensembles rather than one universally described "PDBbind 2019" model. Record the selected model or ensemble and validate it with known co-crystal redocking and, when relevant, cross-docking controls.
