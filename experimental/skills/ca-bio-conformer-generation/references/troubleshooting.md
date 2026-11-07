## ML-Based Conformer Generation and Search

Research methods such as GeoMol generate molecular conformations directly, while TorsionNet uses reinforcement learning to search torsional conformational space. They are distinct approaches and should not be presented as interchangeable drop-in generators.

```python
# Pseudo-code for GeoMol-style ML conformer generation
# (Requires pre-trained model + dependencies)
# from geomol import generate_conformers
# conformers = generate_conformers(smiles, n_conformers=10)
```

**Trade-off:** Performance and coverage depend on the released model, training distribution, conformer definition, and benchmark. Verify the maintained implementation and benchmark it against ETKDGv3 or CREST on the intended chemistry before using it operationally; do not infer broad macrocycle or organometallic coverage from drug-like benchmarks.

## Per-Tool Failure Modes

### ETKDGv3 -- failed embedding

**Trigger:** Macrocycle, highly constrained polycyclic, or sterically crowded molecule.

**Mechanism:** Distance geometry cannot find a consistent 3D structure within the configured embedding iterations.

**Symptom:** `EmbedMolecule` returns -1; `EmbedMultipleConfs` returns empty list.

**Fix:** Set `useRandomCoords=True`, increase `maxIterations`; for macrocycles, set `useMacrocycleTorsions=True`. As fallback, use CREST.

### MMFF94 -- parameter missing

**Trigger:** Molecule contains element not parameterized (transition metals, certain S+ species).

**Mechanism:** MMFF94 only covers H, C, N, O, F, Si, P, S, Cl, Br, I + select cations.

**Symptom:** `MMFFGetMoleculeProperties` returns None; optimization silently no-ops.

**Fix:** Fall back to UFF; or for metals, use GFN2-xTB.

### Conformer ensemble too small

**Trigger:** `n_conf=10` for a flexible molecule (>5 rotatable bonds).

**Mechanism:** A fixed small ensemble can miss relevant minima for flexible molecules.

**Symptom:** New conformers continue to change cluster populations or the downstream result.

**Fix:** As a repository starting heuristic, use n_conf = max(10, 5 * NumRotatableBonds + 10), then increase sampling until the ensemble is stable for the downstream metric.

### Single-conformer 3D descriptor

**Trigger:** Calculating 3D descriptors from a single conformer.

**Mechanism:** Some 3D descriptors vary materially across conformers.

**Symptom:** Same molecule produces different 3D descriptors on rerun.

**Fix:** Converge the descriptor over an ensemble and report the chosen summary. Use Boltzmann weighting only with a justified population model.

### CREST -- timeout on flexible molecule

**Trigger:** Cyclosporin or large peptide.

**Mechanism:** CREST metadynamics scales poorly with rotational complexity.

**Symptom:** Hours of CPU time per molecule; incomplete sampling.

**Fix:** Use `--gfn-ff` for an exploratory lower-cost run or shorten metadynamics with the documented `--mdlen`/`--len` option. `--noopt` disables input pre-optimization; it does not skip metadynamics.

### GFN2-xTB conformer reordering

**Trigger:** Comparing conformer energies between GFN2-xTB and DFT.

**Mechanism:** GFN2-xTB is parameterized for energies; relative conformer ordering can differ from DFT by 1-2 kcal/mol.

**Symptom:** "Wrong" conformer reported as global minimum vs DFT reference.

**Fix:** For high-stakes work, re-rank top GFN2-xTB conformers with DFT single-points (e.g., r2SCAN-3c).

## Reconciliation: ETKDGv3 vs CREST

| Use case | ETKDGv3 | CREST |
|----------|---------|-------|
| Drug-like organic molecule | Efficient starting point | Higher-cost comparison when convergence fails |
| Highly flexible molecule | Increase and convergence-test sampling | Useful alternative sampling strategy |
| Macrocycle or peptide | Try macrocycle-aware settings and validate | Often useful, but not automatically sufficient |
| Population-weighted descriptors | Requires justified energy/population model | Higher-level energy still requires thermodynamic validation |
| FEP input | Useful for initial coordinates | Does not replace bound-pose and MD preparation |

For ETKDGv3 ensembles, compare a representative subset with an orthogonal method or experimental conformers and judge adequacy using the downstream metric; no single cross-method RMSD establishes completeness.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `EmbedMolecule` returns -1 | Embed failed | Set `useRandomCoords=True`; raise `params.maxIterations` |
| MMFFOptimize no-op | MMFF parameters missing | Use UFF fallback |
| All conformers identical | Stiff molecule | OK; molecule is rigid |
| Conformers physically wrong | Stereochemistry lost | Re-add explicit stereo before embedding |
| 3D descriptors differ per run | Random seed not set | `params.randomSeed = 42` |
| CREST out-of-memory | Search/ensemble too large | Reduce `--T`, shorten documented sampling, or lower `--ewin` to retain fewer structures |
| Macrocycle ring inverted | Default torsion preferences wrong | Set `useMacrocycleTorsions=True` |
| AddHs not called | Implicit H not embedded | `mol = Chem.AddHs(mol)` before EmbedMolecule |
