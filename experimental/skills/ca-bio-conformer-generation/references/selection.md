## Conformer Method Taxonomy

| Method | Cost / mol | Quality | Use case | Fails when |
|--------|-----------|---------|----------|------------|
| ETKDGv3 + MMFF94 | Benchmark on actual molecules/hardware | Useful for many drug-like organics | Initial docking/descriptors | Difficult macrocycles, peptides, unsupported chemistry |
| ETKDGv3 + UFF | Fast | Different parameter coverage from MMFF94 | Fallback only after checking UFF parameters | Unsupported atom types; coordination chemistry |
| Omega (OpenEye) | Benchmark licensed workflow | Commercial conformer generator | Commercial pipelines | License cost and configured limits |
| Confab (Open Babel) | Benchmark on intended chemistry | Systematic torsion search | Alternative enumeration | Combinatorial growth and force-field dependence |
| RDKit ETKDGv3 + macrocycle preferences | Molecule-dependent | Macrocycle-aware embedding | Macrocyclic starting ensembles | Coverage remains molecule-dependent |
| CREST + GFN2-xTB | Molecule/settings-dependent | Semiempirical conformational sampling | Difficult flexible molecules | Computational cost; special chemistry |
| CREST + GFN-FF | Lower cost than GFN2-xTB | Force-field-level sampling | Exploratory sampling | Validate coverage and ordering for the chemistry |
| GeoMol (Ganea 2021) | Hardware/model-dependent | Learned conformer generation | Large-library research workflow | Training distribution and released-model coverage |
| TorsionNet (Gogineni 2020) | Hardware/model-dependent | Learned torsional search | Research workflow | Training distribution and implementation availability |
| MD sampling (OpenMM) | System/protocol-dependent | Dynamic sampling | Free energy, induced fit | Computational cost and convergence |

**Decision:** Start drug-like organic molecules with **ETKDGv3** and a parameter-checked MMFF94/MMFF94s optimization. Escalate difficult macrocycles, peptides, or highly flexible molecules to a validated CREST workflow when downstream convergence is inadequate. Benchmark ML generators on the intended chemistry before using them at scale.

## Decision Tree by Scenario

| Scenario | Starting method | Sampling and filtering decision |
|----------|-----------------|---------------------------------|
| Single initial 3D structure | ETKDGv3 + checked force field | Confirm embedding and minimization; downstream relaxation may still be required |
| Multi-conformer docking | ETKDGv3 ensemble | Increase sampling until pose recovery or enrichment is stable |
| 3D descriptors / pharmacophores | ETKDGv3 ensemble | Converge the reported statistic; justify energy/RMSD filters |
| Macrocycle / peptide | Macrocycle-aware ETKDG, then CREST if needed | Compare coverage against known conformers or downstream convergence |
| FEP input | Bound-pose-informed preparation and MD | Do not select solely by isolated-molecule conformer energy |
| Shape search | Query- and library-specific ensemble | Converge retrieval performance on a reference set |
