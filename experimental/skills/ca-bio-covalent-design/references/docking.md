## Covalent Docking Tools

| Tool | Approach | Strength | Fails when |
|------|----------|----------|------------|
| DOCKovalent (London et al 2014 Nat Chem Biol 10:1066) | Constraint-based DOCK | Free, well-validated | Browser-based; small library |
| GOLD covalent (CCDC) | GOLD with covalent constraint | Commercial; selectivity | License cost |
| AutoDock 4 covalent | AD4 with covalent bond | Open source | Slower than Vina |
| CovDock (Schrödinger) | Glide-based + covalent | Commercial two-stage covalent docking workflow | License cost |
| MOE covalent | Triposite Discovery | Commercial | License cost |
| HCovDock (Wu Q, Huang S-Y 2023 Briefings Bioinform 24:bbac559) | Hierarchical fragment + covalent | Open; supports many residues | Newer, less validated |
| ICM-Pro covalent | Active site grid + covalent | Commercial; metal centers | License cost |

For open-source covalent docking, **HCovDock** (2023) is the modern alternative; **DOCKovalent** is the longstanding standard.

**AutoDock 4 covalent, verified runnable (checked 2026-09-21):** the "flexible side chain" method (Bianco G, Forli S, Goodsell DS, Olson AJ, *Protein Sci.* 2016, DOI 10.1002/pro.2733) models the covalent ligand as a flexible receptor residue and docks it with unmodified AutoDock4/AutoGrid4 (`autogrid4 -p <target>.gpf -l <target>.glg`, then `autodock4 -p <ligand>_<target>.dpf -l <ligand>_<target>.dlg`) against `unbound_energy 0.0` in the DPF. Scripps ships AutoDock 4.2.6 Windows binaries (autodock.scripps.edu/download-autodock4) and a covalent-docking tutorial with scripts and a worked PDB 3UPO example (autodock.scripps.edu/resources/covalent-docking); re-running that example's pre-built inputs here reproduced its reference result (best pose -10.73 vs. the shipped -10.74 kcal/mol). Generating the GPF/DPF and flexible-receptor PDBQT for a new target needs MGLTools/AutoDockTools (`prepare_receptor4.py`, `prepare_flexreceptor4.py`, `prepare_gpf4.py`, `prepare_dpf4.py`). The official Scripps Windows installer is an old self-extractor that will not unpack non-interactively, but **bioconda ships a working Linux build** (`bioconda::mgltools` 1.5.7): `micromamba create -n mgltools -c bioconda -c conda-forge mgltools=1.5.7`, run under WSL. **Verified 2026-09-22**: `pythonsh $(which prepare_receptor4.py) -r receptor.pdb -o receptor.pdbqt -A hydrogens` ran clean on the 3PTB fixture (exit 0, 2015-atom PDBQT with Gasteiger charges assigned) — this is a real turnkey pipeline, not engine-only.
