# Per-tool failure modes and Vina/GNINA reconciliation (moved from SKILL.md)

## Per-Tool Failure Modes

### Vina -- cross-dock failure

**Trigger:** Receptor structure not the holo (co-crystal with ligand from another binder).

**Mechanism:** Cross-docking introduces receptor-conformation mismatch, so pose recovery can be substantially worse than self-docking; the size of the decrease is benchmark- and target-dependent.

**Symptom:** Top-ranked pose makes no geometric sense; key contacts missing.

**Fix:** GNINA CNN scoring or ensemble docking. For genuine apo, predict holo with AlphaFold3 / Boltz-1 then dock.

### GNINA CNN -- novel chemotype out-of-distribution

**Trigger:** Ligand chemotype not in PDBbind training.

**Mechanism:** CNN scoring overfits to PDBbind chemotypes; novel macrocycle / peptide / PROTAC scores poorly.

**Symptom:** Affinity prediction far worse than Vina alone.

**Fix:** Use `--cnn_scoring rescore` (sampling still by Vina) rather than CNN sampling. Validate against co-crystal of close analog.

### Box too small

**Trigger:** Binding box defined tightly around small ligand reference.

**Mechanism:** Vina explores only within the box; large analogs cannot fit.

**Symptom:** Many ligands report "no valid pose"; chemotype-biased hits.

**Fix:** Derive the box from the reference ligand or known pocket and add enough explicit padding for the largest intended ligands to translate and rotate. Then verify containment and redocking/search convergence on controls. There is no universal padding value or 25 A cube that fits every ligand series.

### Multi-pocket protein -- wrong site

**Trigger:** Protein has multiple binding sites (orthosteric + allosteric).

**Mechanism:** P2Rank or AutoBox picks the most "drugable" pocket; not always the desired one.

**Symptom:** Hits dock in wrong pocket; SAR confusing.

**Fix:** Verify pocket from co-crystal data; explicitly set `center_x/y/z` from known ligand centroid.

### DiffDock-L -- PoseBusters invalid

**Trigger:** Default DiffDock-L output for any receptor.

**Mechanism:** Diffusion-generated poses are not guaranteed to satisfy every bond-geometry, stereochemistry, and intermolecular-clash check; failure rates vary by method and benchmark.

**Symptom:** Poses look reasonable but fail PoseBusters checks.

**Fix:** Filter to PB-valid (PoseBusters); rescore with GNINA. See `chemoinformatics/pose-validation`.

### Wrong ionization state

**Trigger:** Ligand or receptor residues protonated incorrectly at pH 7.4.

**Mechanism:** Aspartate/glutamate/histidine protonation depends on local environment; default protonation may be wrong.

**Symptom:** Salt bridges missing; poses misranked.

**Fix:** Run PROPKA on the receptor to estimate residue pKas; for catalytic histidines, manually inspect protonation and tautomer state in the local environment.

## Reconciliation: Vina vs GNINA Disagreement

| Vina top pose | GNINA top pose | Action |
|---------------|----------------|--------|
| Same pose, similar score | Same pose, similar score | Treat agreement as supporting evidence; still run physical-validity checks |
| Vina top pose ≠ GNINA top pose | Same pocket, different orientation | Retain both and compare against target-relevant controls or interaction evidence |
| Vina excellent, GNINA mediocre | Different pose, very different score | Inspect both poses; do not infer which method is correct from score disagreement alone |
| Both poor scores | Many ligands score similarly poor | Wrong pocket / protein conformation; reconsider receptor |
