## PoseBusters Test Suite

PoseBusters runs ~20 individual checks grouped into:

The thresholds below are the benchmark criteria reported by Buttenschoen et al. (2024). Installed PoseBusters defaults may differ by version and configuration, so record the package version and resolved configuration.

| Check group | What it tests | 2024 benchmark criterion |
|-------------|---------------|-----------|
| Sanity | Ligand chemical sanity | RDKit sanitization passes |
| Bond lengths | Bond lengths within reference | 0.75–1.25 times RDKit distance-geometry bounds |
| Bond angles | 1–3 distances within reference | 0.75–1.25 times RDKit distance-geometry bounds |
| Internal steric | No intra-ligand clash | Pair distance > 0.70 times the RDKit lower bound |
| Aromatic ring planarity | Aromatic rings planar | Maximum deviation from fitted plane <= 0.25 Å |
| Double-bond stereo | Z/E preserved | Match input SMILES |
| Internal energy | Energy relative to generated conformers | UFF energy ratio <= 100 versus the mean of 50 generated, relaxed conformers |
| Volume overlap | vdW overlap with protein | < 7.5% of ligand vdW volume |
| Minimum distance | No severe protein-ligand clash | Distance >= 0.75 times the sum of vdW radii |
| Chirality | R/S preserved from input | Match input SMILES |

Double-bond stereo, chirality, molecular formula, and molecular bond identity are reference checks: they only run when `mol_true` is supplied (`config='redock'`). Without a reference (`config='dock'`), those checks -- and RMSD -- are absent from the results, not merely unreported (verified against PoseBusters 0.6.5: `dock` returns 22 boolean columns, `redock` returns 28, the difference being `rmsd_<=_2å`, `mol_true_loaded`, `molecular_formula`, `molecular_bonds`, `double_bond_stereochemistry`, and `tetrahedral_chirality`).

A pose passing ALL tests is "PB-valid". Combined PB-valid + RMSD <= 2 Å is the modern criterion.

## When to Apply PoseBusters

| Workflow | PoseBusters use | Action |
|----------|-----------------|--------|
| Self-docking (validating method) | Required | Compare PB-valid + RMSD <= 2A |
| Cross-docking | Required | PB-valid + RMSD <= 2A; account for protein flexibility |
| Virtual screening top hits | Required | Filter to PB-valid before MM/GBSA / FEP |
| AI docking (DiffDock, etc.) | Required for a fair benchmark | Report the dataset-specific PB-valid and combined success rates |
| Generated ligand poses | Recommended | Measure chemical and geometric validity rather than assuming it |
| Boltz-2 / AlphaFold3 ligand poses | Recommended | Benchmark validity on the relevant complexes; do not infer a failure frequency from DiffDock |
| Production FEP setup | Required | Inspect pose validity and ligand strain before system preparation |

## PoseBusters Usage

```python
from posebusters import PoseBusters

bust = PoseBusters(config='redock')

results = bust.bust(
    mol_pred='predicted.sdf',
    mol_true='reference.sdf',
    mol_cond='receptor.pdb',
)
```

Common configurations and their included checks are:

| Config | Includes | When to use |
|--------|----------|-------------|
| `redock` | All checks + RMSD vs reference + protein vdW overlap | Self-docking benchmarks, retrospective validation |
| `dock` | All checks that need no reference; drops RMSD and the reference-dependent checks listed under the check table | Blind docking, prospective virtual screening |
| `mol` | Intra-ligand only (sanitization, bond lengths/angles, internal clash, ring and double-bond flatness, energy); no stereo or chirality checks, which need `mol_true` | Conformer QC; no protein context |

PoseBusters also ships additional and faster configurations in some releases. Treat the table as a workflow guide, not an exhaustive registry, and inspect the configurations available in the installed version.

Output: a DataFrame with one row per pose, metadata columns, and boolean pass/fail columns for the checks enabled by the selected configuration. Reference-dependent fields such as RMSD and the exact check-column names vary by configuration and version; inspect `results.columns` rather than relying on a fixed exhaustive list.

## Python Library API

**Goal:** Programmatically validate a docked-pose SDF against a receptor PDB and produce a PB-valid filter.

**Approach:** Instantiate `PoseBusters(config='dock')`, call `bust()` on the SDF + PDB pair, and AND-aggregate all boolean check columns into a single `pb_valid` flag.

Run `run_posebusters(pred_sdf, receptor_pdb)` from `scripts/validate_poses.py` (it uses `dock`, or `redock` when a reference SDF is given, and adds the `pb_valid` column). Check columns are the boolean ones excluding `rmsd*`.
