---
name: bio-virtual-screening
description: Performs structure-based virtual screening using AutoDock Vina, SMINA, GNINA (CNN scoring), and DiffDock-L hybrid workflows with explicit choice rules across rigid vs flexible docking, cross-docking vs self-docking, binding-site detection (P2Rank, fpocket), receptor preparation (PDB2PQR, PROPKA), ligand preparation (meeko, OpenBabel), and ultralarge-library screening (ZINC22, Enamine REAL). Use when screening chemical libraries against a protein target to find candidate binders, ranking docking poses, or selecting a docking workflow for a specific scenario.
tool_type: python
primary_tool: AutoDock Vina
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: AutoDock Vina 1.2.5+, SMINA 2020-12+, GNINA 1.1+ for `rescore` (GNINA 1.3+ for the six-mode interface documented below), RDKit 2024.09+, meeko 0.5+, P2Rank 2.4+, ProDy 2.4+, pdb2pqr 3.6+. Receptor prep verified end-to-end (pdb2pqr 3.7.1, meeko 0.8.0) against PDB 3PTB. Vina 1.1.2 vs 1.2 may give different poses for the same input -- check the installed major version before comparing runs.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `vina --version`; `gnina --version`; `smina --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

**meeko's CLI entry point has no `.py` suffix.** `pip install meeko` (0.5+) registers a
console-script named `mk_prepare_receptor` and `mk_prepare_ligand` -- not
`mk_prepare_receptor.py`. The `.py`-suffixed form raises `FileNotFoundError` once
installed via pip; use the bare name, as in the code below.

**On Windows, `pip install vina` has no wheel** (`ValueError: Boost library location was
not found!` at build time). Use the Vina CLI via `subprocess` instead of `from vina import
Vina` -- see the CLI-fallback comment in "Vina Docking (Single Ligand)" below. The Python-API
branch itself is verified, not just doc-checked: conda-forge ships prebuilt `vina` 1.2.7 for
linux-64 (`conda install -c conda-forge vina`), which runs under WSL on a Windows box. Docking
benzamidine into PDB 3PTB there with `scripts/dock_single.py`'s `from vina import Vina` path and
with `examples/virtual_screen.py`'s `virtual_screen()` gave the same top pose (-5.978 kcal/mol)
as the Vina CLI run on the same fixture, and `write_poses`/`v.energies()` round-tripped a valid
multi-model PDBQT. No Windows build of the PyPI `vina` wheel exists, so the CLI fallback stays
the Windows-native path -- but the API branch is confirmed correct, not merely plausible.

# Virtual Screening

Screen chemical libraries against protein targets via molecular docking. Vina is the de-facto default, SMINA adds flexibility (Vinardo scoring, custom scoring), and GNINA adds CNN-based pose scoring (Top-1 redock 58%->73% over Vina, cross-dock 27%->37%). Deep-learning docking (DiffDock-L, EquiBind, NeuralPLexer) competes in pose accuracy, but physical validity is method- and dataset-dependent; the workflow therefore combines ML pose sampling with classical scoring and explicit geometry checks. For ultralarge libraries (>1M), library preparation, hierarchical filtering, and HPC orchestration become the limiting steps.

For pose physical-validity QC, see `chemoinformatics/pose-validation`. For ML-driven docking + rescoring, see `chemoinformatics/ml-docking-rescoring`. For covalent docking, see `chemoinformatics/covalent-design`. For affinity calculations (FEP), see `chemoinformatics/free-energy-calculations`.

**Handoff caveat:** Open Babel's PDBQT -> SDF conversion (`obabel out.pdbqt -O pose.sdf`) loses formal bond order and charge for charged ligands (PDBQT encodes neither), so RDKit sanitization fails and PoseBusters passes only 3/12 checks even when the pose is spatially valid. Rebuild the pose from the SMILES that meeko writes into the PDBQT instead (Vina keeps the `REMARK SMILES` lines in its output). Verified on benzamidine docked to trypsin (3PTB): `obabel` route 3/12, meeko route 12/12 on all 5 poses (meeko 0.8.0, RDKit 2026.03.6, PoseBusters 0.6.5, Vina 1.2.7):

```python
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem

pm = PDBQTMolecule.from_file('out.pdbqt', skip_typing=True)   # Vina output, all poses
mol = RDKitMolCreate.from_pdbqt_mol(pm)[0]                    # one Mol, one conformer per pose
with Chem.SDWriter('poses.sdf') as w:
    for cid in range(mol.GetNumConformers()):
        w.write(mol, confId=cid)
```

Then `bust poses.sdf --outfmt short`. This needs the ligand PDBQT to come from meeko (`prepare_ligand`), which writes the `REMARK SMILES` lines; a PDBQT from another writer has none, so carry the original RDKit `Mol` alongside it instead.

## Docking Tool Taxonomy

| Tool | Scoring | Speed (sec/lig) | Best at | Fails when |
|------|---------|-----------------|---------|------------|
| AutoDock Vina 1.2 | Vina (empirical) | Hardware- and settings-dependent | Open, well-characterized baseline | Cross-dock; cryptic pockets; metal centers |
| SMINA | Vina + flexible + custom | Hardware- and settings-dependent | Custom scoring; flexible side chains | Same Vina-scoring caveats |
| Vinardo | Modified Vina scoring | Hardware- and settings-dependent | Alternative empirical score | Validate on target-relevant controls |
| GNINA 1.1 | CNN or Vina scoring | GPU- and settings-dependent | CNN-assisted pose ranking | Validate transfer to the target and chemotype |
| AutoDock 4 | AD4 + grid maps | Hardware- and settings-dependent | Legacy reference | More setup than Vina |
| DOCK 6/7 | DOCK + Amber | Hardware- and settings-dependent | UCSF DOCK ecosystem | Steep learning curve |
| Glide (Schrodinger) | GlideScore | License and hardware-dependent | Commercial docking workflow | License cost |
| GOLD (CCDC) | GOLDScore / ChemScore | License and hardware-dependent | Commercial workflow; metal options | License cost |
| FlexX (BioSolveIT) | FlexX | License and hardware-dependent | Fragment-based placement | License cost |
| rDock | rDock | Hardware- and settings-dependent | Open-source alternative | Validate maintenance and target fit |
| DiffDock-L | Diffusion-generative | GPU- and settings-dependent | Pose sampling for cross-docking | Validate geometry with PoseBusters; see ml-docking-rescoring |
| EquiBind | Equivariant NN | GPU- and settings-dependent | Single-shot pose generation | Requires independent geometry and ranking checks |
| Boltz-2 + GNINA rescore | Foundation model + CNN | GPU- and settings-dependent | Experimental multi-model workflow | Benchmark each evidence stream independently |

**Decision:** Use Vina as an open baseline and consider GNINA CNN rescoring when target-relevant redocking or cross-docking controls support it. For large libraries, calibrate a hierarchical Vina -> GNINA -> higher-cost follow-up workflow on measured enrichment, throughput, and retained chemotype diversity.

## Decision Tree by Scenario

| Scenario | Recommended workflow |
|----------|---------------------|
| Self-dock against known ligand pocket | GNINA `gnina --cnn_scoring rescore` (see `references/gnina.md`) |
| Cross-dock to apo or related-target structure | DiffDock-L pose + GNINA rescore + PoseBusters (failure modes: `references/failure-modes.md`) |
| Ultralarge library (10M+) | Calibrated hierarchical screen: property/alert triage -> Vina -> measured top fraction to GNINA -> higher-cost follow-up (see `references/ultralarge-screening.md`) |
| Cryptic pocket / induced fit | Receptor-ensemble docking and, where appropriate, a separately validated complex-prediction model |
| Allosteric / undefined site | P2Rank for pocket detection -> ensemble dock all pockets |
| Metal-coordinated ligand | GOLD (commercial) or manually parameterize Vina metal scoring |
| Covalent inhibitor | See `chemoinformatics/covalent-design`: DOCKovalent, HCovDock |
| Fragment screen (<300 Da) | rDock or constrained Vina with seed atoms |
| Hit-to-lead refinement | Use co-crystal structure if available; MD-relaxed receptor; FEP for affinity |

## Receptor Preparation

**Goal:** Convert a protein PDB into a docking-ready format with correct protonation, missing atoms, and removed waters.

**Approach:** Decide which ligands, cofactors, metals, and structural waters to retain -> fill missing heavy atoms with a structure-repair tool such as PDBFixer -> use PROPKA/PDB2PQR plus manual review to assign pH-dependent protonation -> assign the charge model required by the docking workflow -> prepare receptor PDBQT with a documented AutoDock-compatible tool.

```bash
python scripts/prepare_receptor.py repaired.pdb receptor.pdbqt --ph 7.4
```

`scripts/prepare_receptor.py` runs `pdb2pqr --ff=AMBER --with-ph=<pH> --pdb-output` then `mk_prepare_receptor --read_pdb ... -p`. Decide which waters, cofactors and metals to retain before calling it. It needs `pdb2pqr` and `mk_prepare_receptor` on `PATH`; it writes the protonated `<base>_pH<pH>.pdb`/`.pqr` beside the input. It hands meeko the protonated PDB rather than the `.pqr` (see the insertion-code pitfall below).

**Common pitfall:** Forgetting to add hydrogens at protein pH (7.4) but using pH 7.0 ligand charges. Hist mistakenly protonated. Use PROPKA + manual review of catalytic residues.

**Common pitfall:** Feeding pdb2pqr's default `.pqr` output straight into `mk_prepare_receptor --read_pqr`. meeko 0.8.0's PQR reader assumes an all-integer residue-number column and raises `ValueError: invalid literal for int() with base 10` on any residue with a PDB insertion code (e.g. `184A`). Any PDB deposition with insertion-code residues triggers it, not one protein family: chymotrypsin-numbered serine proteases (trypsin 3PTB: 184A/188A/221A; elastase 1EAI) are common examples and hit it on real structures, not just edge cases. Use `pdb2pqr --pdb-output` and `mk_prepare_receptor --read_pdb` as above; verified on PDB 3PTB (trypsin, insertion-code residues 184A/188A/221A).

## Ligand Preparation

**Goal:** Generate a 3D, docking-ready ligand file from SMILES with appropriate protonation and conformation.

**Approach:** Supply a documented protomer/tautomer state generated by an appropriate pKa/protomer workflow -> parse it with RDKit -> embed 3D with ETKDGv3 -> minimize with MMFF94 -> write PDBQT with Meeko. `MolFromSmiles` parses the supplied state and `Uncharger` neutralizes formal charges; neither predicts protonation at pH 7.4.

```python
import sys; sys.path.insert(0, 'examples')      # this skill's examples/ directory
from virtual_screen import prepare_ligand
prepare_ligand('NC(=[NH2+])c1ccccc1', 'ligand.pdbqt')   # supplied protomer SMILES -> ETKDGv3 -> MMFF94 -> meeko PDBQT
```

`prepare_ligand` in `examples/virtual_screen.py` raises on invalid SMILES, embedding failure, missing MMFF94 parameters, MMFF94 non-convergence, or a meeko export error (meeko 0.5+ API: `MoleculePreparation().prepare(mol)` returns setups; `PDBQTWriterLegacy.write_string(setups[0])` materializes the PDBQT).

`meeko` (AutoDock developers' tool) handles torsion tree creation, rotamer flagging, and PDBQT writing -- preferred over Open Babel's PDBQT writer. Note: meeko 0.5+ separated the writer (`PDBQTWriterLegacy`) from `MoleculePreparation`; older code using `prep.write_pdbqt_file()` is deprecated.

## Binding Site Detection

When the binding pocket is not known (apo target, novel allosteric site):

| Tool | Approach | Output |
|------|----------|--------|
| P2Rank (Krivak 2018) | ML on protein surface descriptors | Ranked pocket list with center coords |
| fpocket (Le Guilloux 2009) | Voronoi tessellation | Pocket descriptor list |
| DoGSiteScorer | Geometric + drugability | Pocket list with score |
| AutoSite (Vina) | Affinity map clustering | Pocket centers |
| AlphaFill | Transplant ligands/cofactors from homologous experimental structures into AlphaFold models | Plausible binding-site components for review |

```bash
prank predict -f receptor.pdb -o pockets/
```

**On Windows, `prank.bat` returns immediately with no output when invoked through a shell
wrapper** (`cmd.exe /c`, or a bash script calling it) -- confirmed on this fixture: it printed
only the `cmd.exe` banner and produced no `pockets/` directory. Call the jar directly instead,
which does run and predict:

```bash
java -cp "<p2rank_dir>/bin/p2rank.jar;<p2rank_dir>/bin/lib/*" cz.siret.prank.program.Main \
     predict -f receptor.pdb -o pockets/
```

Verified on PDB 3PTB: 4 pockets in ~4 s, `pockets/rec.pdb_predictions.csv` with `pocket1` centered
at (-1.52, 14.47, 17.47) -- the same box center used throughout this Skill's Vina examples, which
sits on the real benzamidine site. P2Rank output `<receptor>_predictions.csv` lists pocket centers
with scores. The highest model score does not identify a pocket as orthosteric or biologically
relevant; verify ranked pockets against co-crystal, mutagenesis, SAR, or other structural evidence.

## Vina Docking (Single Ligand)

```bash
python scripts/dock_single.py receptor.pdbqt ligand.pdbqt --center X Y Z --size X Y Z \
       --exhaustiveness 8 --n-poses 10 --seed 42 --out poses.pdbqt
```

`scripts/dock_single.py` uses the Vina Python API (Vina 1.2+) when `from vina import Vina` imports, and otherwise runs the Vina CLI (`--vina-exe` if `vina` is not on `PATH`; Windows has no `vina` wheel, see Version Compatibility). Both paths pass the same seed. It writes the poses to `--out` and prints one affinity per mode with `affinity >= 0` modes dropped (Vina occasionally emits a physically nonsensical positive-energy mode, e.g. +68 kcal/mol; see "Sanity-filter reported poses").

**Exhaustiveness:** `8` is the Vina default. Increasing it increases search effort, but runtime and pose recovery depend on hardware, ligand flexibility, box size, and software version. Benchmark settings such as 8, 16, 32, and 64 on target-relevant controls instead of assigning universal timing or quality labels.

**Seed and reproducibility:** `scripts/dock_single.py` and `examples/virtual_screen.py` accept a `seed` (CLI: `--seed`). Top-1 affinity is empirically stable run-to-run without a fixed seed, but poses ranked 2+ reorder between runs on the same input. Set and record a seed (default `42` above) whenever the top-N poses -- not only the single best -- will be reported or compared.

**Sanity-filter reported poses:** Vina's raw mode list can include a physically nonsensical outlier (a positive-energy mode was observed among 9 returned modes in testing). Filter or flag `affinity >= 0` poses before reporting or ranking; do not assume every mode Vina returns is a plausible binder.

Vina's `rmsd_lb` and `rmsd_ub` are lower and upper heavy-atom RMSD bounds between a reported mode and the best-scoring mode; the bounds differ in how symmetry-equivalent atoms are handled. They are not pose-versus-experimental-reference RMSDs. Use an external symmetry-aware RMSD to a reference pose for accuracy QC.

## Reference Files

| File | Read when |
|------|-----------|
| `references/gnina.md` | Running GNINA: CLI, `--cnn_scoring` modes, choosing a CNN model or ensemble |
| `references/ultralarge-screening.md` | Screening more than ~100k compounds: hierarchical Vina -> GNINA pipeline skeleton, ZINC22/Enamine REAL access, stage-by-stage heuristics |
| `references/failure-modes.md` | A docking result looks wrong (cross-dock, box too small, wrong pocket, wrong ionization, DiffDock-L invalid poses, GNINA out-of-distribution), or Vina and GNINA disagree |

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

## References

- Trott & Olson, *J. Comput. Chem.* 31:455-461 (2010) -- AutoDock Vina. https://doi.org/10.1002/jcc.21334
- Eberhardt et al., *J. Chem. Inf. Model.* 61:3891-3898 (2021) -- Vina 1.2 features. https://doi.org/10.1021/acs.jcim.1c00203
- AutoDock Vina, official manual -- result-field and CLI semantics. https://vina.scripps.edu/manual/
- AutoDock Vina, official zinc-metalloprotein tutorial -- AutoDock4Zn maps through the Vina executable. https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html
- Quiroga & Villarreal, *PLoS ONE* 11:e0155183 (2016) -- Vinardo scoring. https://doi.org/10.1371/journal.pone.0155183
- McNutt et al., *J. Cheminformatics* 13:43 (2021) -- GNINA 1.0 CNN docking. https://doi.org/10.1186/s13321-021-00522-2
- GNINA, official repository -- current CLI modes and named CNN ensembles. https://github.com/gnina/gnina
- Buttenschoen et al., *Chem. Sci.* 15:3130-3139 (2024) -- PoseBusters benchmark. https://doi.org/10.1039/D3SC04185A
- Lyu et al., *Nature* 566:224-229 (2019) -- ultralarge virtual-screening proof of concept. https://doi.org/10.1038/s41586-019-0917-9
- Krivak & Hoksza, *J. Cheminformatics* 10:39 (2018) -- P2Rank. https://doi.org/10.1186/s13321-018-0285-8
- Forli et al., *Nat. Protoc.* 11:905-919 (2016) -- AutoDock suite and AutoDockTools. https://doi.org/10.1038/nprot.2016.051
- Le Guilloux, Schmidtke & Tuffery, *BMC Bioinformatics* 10:168 (2009) -- fpocket. https://doi.org/10.1186/1471-2105-10-168
- Meeko, official documentation -- ligand/receptor PDBQT preparation and export interfaces. https://meeko.readthedocs.io/
- Dolinsky et al., *Nucleic Acids Res.* 35:W522-W525 (2007) -- PDB2PQR. https://doi.org/10.1093/nar/gkm276
- Olsson et al., *J. Chem. Theory Comput.* 7:525-537 (2011) -- PROPKA 3. https://doi.org/10.1021/ct100578z
- PDBFixer, official repository -- missing-residue/atom repair interface. https://github.com/openmm/pdbfixer
- Corso et al., *ICLR* (2024) -- DiffDock-L. https://proceedings.iclr.cc/paper_files/paper/2024/file/db334db287337b2a365120b524300ef3-Paper-Conference.pdf
- Stärk et al., *ICML* (2022) -- EquiBind. https://proceedings.mlr.press/v162/stark22b.html
- Qiao et al., *Nat. Mach. Intell.* 6:195-208 (2024) -- NeuralPLexer. https://doi.org/10.1038/s42256-024-00792-z
- Passaro et al., bioRxiv (2025) -- Boltz-2. https://doi.org/10.1101/2025.06.14.659707
- Abramson et al., *Nature* 630:493-500 (2024) -- AlphaFold 3. https://doi.org/10.1038/s41586-024-07487-w
- ZINC22, official resource. https://zinc22.docking.org/
- Enamine REAL, official resource. https://enamine.net/compound-collections/real-compounds
- ChEMBL, official versioned database. https://www.ebi.ac.uk/chembl/

## Related Skills

- chemoinformatics/molecular-io - Parse ligands
- chemoinformatics/conformer-generation - Generate 3D for ligand prep
- chemoinformatics/molecular-standardization - Canonicalize before docking
- chemoinformatics/pose-validation - PoseBusters physical-validity QC
- chemoinformatics/ml-docking-rescoring - DiffDock-L + GNINA hybrid
- chemoinformatics/covalent-design - Covalent docking
- chemoinformatics/free-energy-calculations - FEP for refined affinity
- chemoinformatics/admet-prediction - Filter library before docking
- structural-biology/structure-io - PDB / mmCIF handling
- structural-biology/modern-structure-prediction - AlphaFold3 / Boltz-1 for apo receptors
