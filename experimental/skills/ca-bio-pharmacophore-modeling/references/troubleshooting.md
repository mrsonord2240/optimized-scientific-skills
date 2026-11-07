## Pharmacophore vs Shape vs 2D Fingerprint

| Method | Captures | Best for |
|--------|----------|----------|
| ECFP4 Tanimoto | Local atom environments | Lead optimization (same series) |
| FCFP4 Tanimoto | Pharmacophore-equivalent atoms | Loose similarity in series |
| Shape similarity (ROCS) | 3D shape volume | Scaffold hopping by shape |
| Pharmacophore | Discrete features in space | Scaffold hopping with feature specificity |
| Combined (Tanimoto + shape) | Multi-objective | Production VS |

Pharmacophore is more *interpretable* than shape: a hit explains why it matched (donor at position X, hydrophobe at position Y).

## Per-Tool Failure Modes

### Ligand-based -- diverse actives confound

**Trigger:** Active set spans multiple scaffolds with different bound conformations.

**Mechanism:** No common pharmacophore exists; algorithm forces non-consensus features.

**Symptom:** Pharmacophore matches no actives in retrospective.

**Fix:** Cluster actives by scaffold first; derive per-cluster pharmacophore. For the coarse feature-family prefilter in `scripts/pharmacophore.py`, `feature_family_prefilter(..., min_shared_fraction=0.8)` also lets a molecule lacking a few non-essential shared feature types through (default 1.0 is strict); it stays a coarse filter, so follow with a distance-constrained 3D match.

### Receptor-based -- apo structure

**Trigger:** Protein in apo form (no bound ligand).

**Mechanism:** Side-chain rotamers differ between apo and holo; "binding site" geometry is wrong.

**Symptom:** Pharmacophore inferred from apo doesn't match holo experimental data.

**Fix:** Use AlphaFold3 / Boltz-1 to predict holo conformation; derive pharmacophore from predicted holo.

### Pharmacophore -- single conformer bias

**Trigger:** Active aligned to its first generated conformer, not bioactive conformer.

**Mechanism:** Crystal structure not available; generated conformer may not be the bound one.

**Symptom:** Pharmacophore inconsistent across runs (different starting conformer chosen).

**Fix:** Use conformer ensemble; align all to common scaffold; choose conformer most consistent with other actives.

### Tolerance too tight

**Trigger:** Default geometric tolerance < 0.5 Å.

**Mechanism:** Real bioactive conformers have flexibility; rigid pharmacophore filters most molecules out.

**Symptom:** Search returns zero hits.

**Fix:** Use tolerance 1.0-1.5 Å for drug-like; up to 2 Å for flexible peptide-like.

### Pharmacophore search misses bioisostere

**Trigger:** Bioisostere replacement (e.g., -COOH replaced by tetrazole).

**Mechanism:** Tetrazole functions as acid bioisostere but RDKit features may not classify identically.

**Symptom:** Known bioisosteric active not found.

**Fix:** Use ChemAxon-style bioisosteric feature equivalence; or pharmacophore feature class expansion (acid generic vs -COOH specific).

### PLIP -- water bridge absent from output

**Trigger:** Bridging water between ligand donor and protein acceptor.

**Mechanism:** PLIP can report water bridges, but the required crystallographic water must be present in the input and satisfy its geometric criteria.

**Symptom:** Pharmacophore missing critical H-bond feature.

**Fix:** Retain relevant crystallographic waters, inspect PLIP water-bridge output, and review borderline geometry manually.

## Reconciliation: Ligand-Based vs Receptor-Based

| Aspect | Ligand-based | Receptor-based |
|--------|--------------|----------------|
| Data needed | Multiple actives with defensible conformers/alignment | A defined pocket, optionally with a co-crystal ligand |
| Main bias | Known active chemotypes, conformer choice, and alignment | Pocket structure, protonation, retained waters, and interaction-detection/modeling rules |
| Hit-set behavior | Depends on feature abstraction and tolerances | Depends on selected pocket interactions, excluded volumes, and tolerances |
| Confidence evidence | Retrospective recovery across held-out actives/inactives | Recovery of known interaction geometry and retrospective or prospective validation |

Choose between ligand- and receptor-based models using the available structural/activity evidence and target-relevant validation. Neither approach is universally more reliable, diverse, or suitable for scaffold hopping.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `Pharm3D.EmbedPharmacophore` fails | Bounds matrix infeasible | Review/loosen justified bounds and, when more attempts are warranted, increase the documented `count` argument; inspect `n_failed` |
| Pharmacophore matches everything | Too few features | Add features; tighten tolerances |
| Pharmacophore matches nothing | Too many features or tight bounds | Reduce feature count; loosen tolerances |
| BaseFeatures.fdef not found | RDKit installation issue | Check `from rdkit.RDPaths import RDDataDir` |
| Pharmacophore-conformer mismatch | Wrong conformer used | Use bioactive conformer from crystal |
| Pharmit search timeout | Library too large | Pre-filter by 2D fingerprint Tanimoto |
| apo2ph4 PML has no useful model | No robust pocket hot spots at selected settings | Recheck pocket definition and documented thresholds; inspect alternative models |
| PLIP `analyze()` raises `ValueError: inchikey is not a recognised Open Babel format` | `openbabel-wheel` (the standard pip build, checked 3.1.1.23 on Windows) ships no InChI writer at all, and PLIP's `Ligand.__init__` unconditionally calls `pybel.write(format='inchikey')` | Use a conda-forge Open Babel build with InChI support, or monkeypatch `pybel.Molecule.write` to no-op for `format='inchikey'` |
