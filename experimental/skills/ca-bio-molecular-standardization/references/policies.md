## Standardization Pipeline Stages

| Stage | RDKit Tool | Operation | Common errors caught |
|-------|-----------|-----------|----------------------|
| 1. Sanitization | `Chem.SanitizeMol` | Kekulize, assign aromaticity, fix valences | Wrong valence on N/O |
| 2. Salt stripping | `rdMolStandardize.FragmentRemover` or `LargestFragmentChooser` | Remove counterions | Cl-, Na+, K+, OH- |
| 3. Mixture choice | `LargestFragmentChooser` | Pick parent fragment | Co-crystals, hydrates |
| 4. Charge neutralization | `Uncharger` | Neutralize while preserving net charge | Permanent charges preserved (quaternary N+) |
| 5. Tautomer canonicalization | `TautomerEnumerator.Canonicalize` | Pick canonical tautomer | Keto/enol; amide/imidate |
| 6. Stereo standardization | `Chem.AssignStereochemistry` | Consistent stereo descriptors | Lost wedges, ambiguous R/S |
| 7. Isotope normalization | Explicitly set selected atom isotope labels to 0 | Remove 13C, 2H labels | Tracer studies; preserve labels when scientifically meaningful |
| 8. Output canonicalization | `Chem.MolToSmiles(canonical=True)` | Canonical SMILES + InChIKey | Round-trip stability |

## Pipeline Reconciliation

| Pipeline | Origin | Tautomer canonicalization | Salt definition | Use case |
|----------|--------|---------------------------|-----------------|----------|
| ChEMBL pipeline | EBI ChEMBL | Not performed by `standardize_mol` or `get_parent_mol` | ChEMBL salt list (extensive) | ChEMBL-compatible registration |
| canSARchem | ICR Cancer Research UK | Canonical tautomer BEFORE parent extraction | Extended salt list | Cancer drug discovery |
| PubChem (OpenEye) | NIH NCBI | OpenEye QUACPAC tautomer | PubChem salt list | Bioassay data, large-scale |
| RDKit rdMolStandardize default | Greg Landrum | RDKit TautomerEnumerator | RDKit default | General purpose, open source |

**Key difference (canSARchem vs ChEMBL):**
- ChEMBL standardizes the representation and extracts a parent, but does not canonicalize tautomers.
- canSARchem canonicalizes the tautomer before parent extraction.

This difference matters when alternate tautomeric inputs must be registered as one parent. Do not describe ChEMBL output as tautomer-canonical unless an explicit tautomer step is added and documented.
