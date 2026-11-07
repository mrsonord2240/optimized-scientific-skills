## Per-Tool Failure Modes

### ChEMBL pipeline -- multi-fragment parent

**Trigger:** Every fragment is on ChEMBL's salt list, including an organic salt such as sodium acetate, sodium benzoate, sodium formate, sodium citrate, potassium acetate, or choline chloride.

**Mechanism:** `get_parent_mol` can return its input unchanged because no fragment qualifies as a surviving parent under its salt-list rules.

**Symptom:** A supposedly standardized ChEMBL parent still has two or more fragments. A carbon-count pre-filter does not catch this: these examples contain carbon.

**Fix:** Check `len(Chem.GetMolFrags(parent)) > 1` immediately after `get_parent_mol`. The safe default is to flag the record for review. Use `LargestFragmentChooser(preferOrganic=True)` only as an explicit fallback, and record that it may discard meaningful co-crystal information.

### ChEMBL pipeline -- inorganic parent

**Trigger:** Molecule is genuinely an inorganic salt (e.g., NaCl, K2SO4).

**Mechanism:** `get_parent_mol` chooses largest organic; falls back to largest fragment for fully inorganic.

**Symptom:** Returns the salt itself (not a drug).

**Fix:** After the fragment-count check, reject/flag a selected parent with no carbon atoms. This is a separate condition from the all-fragments-are-salts case.

### Uncharger -- charge-state policy mismatch

**Trigger:** A molecule combines a non-removable charge, such as quaternary ammonium, with other neutralizable sites, or the desired physiological ionization state differs from a structure-normalization rule.

**Mechanism:** `Uncharger` adds or removes hydrogens from neutralizable acids and bases. It cannot remove a permanent charge that has no corresponding hydrogen edit; by default it may preserve an opposite neutralizable charge when a non-removable charge is present so that the total charge remains balanced. `force=True` instead neutralizes all sites that can be neutralized even if the remaining permanent charge leaves a nonzero total charge.

**Symptom:** The permanent charge remains, but other sites or the total charge differ from the protonation state intended for docking or modeling.

**Fix:** Choose `force` according to the documented total-charge policy, keep `force=False` when balanced countercharges should be preserved, and inspect/prepare physiological protonation states separately.

### Tautomer enumerator -- combinatorial explosion

**Trigger:** Molecule with many tautomerizable groups (polyhydroxylated heterocycle).

**Mechanism:** `TautomerEnumerator.Enumerate` generates all possible tautomers; can produce thousands.

**Symptom:** OOM or hour-long compute on single molecule.

**Fix:** Use `Canonicalize` when only the configured canonical representation is needed. Before `Enumerate`, call `enumerator.SetMaxTransforms(limit)` (and, when appropriate, `SetMaxTautomers(limit)`) to cap the search.

### Legacy MolVS -- import or compatibility failure

**Trigger:** Code still using legacy `from molvs import Standardizer`.

**Mechanism:** The standalone MolVS package is legacy and may not support current Python/RDKit versions. RDKit's maintained `rdMolStandardize` module remains available.

**Symptom:** ImportError or AttributeError on newer RDKit.

**Fix:** Migrate deliberately to `from rdkit.Chem.MolStandardize import rdMolStandardize`; compare outputs because RDKit functions are not drop-in aliases for every MolVS workflow.

### Round-trip InChIKey mismatch

**Trigger:** Records were processed with different standardization settings or entered in different salt, charge, isotope, stereo, or tautomer forms.

**Mechanism:** The pipelines did not apply the same explicitly versioned transformations before identity generation.

**Symptom:** Apparently equivalent records produce different InChIKeys, or an expected database join fails.

**Fix:** Record and apply the same toolkit version, standardization stages, tautomer policy, and InChI options to both datasets; compare full standardized structures when results still differ.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| ImportError from standalone `molvs` | Legacy package incompatible with current environment | Use maintained `rdkit.Chem.MolStandardize.rdMolStandardize` APIs and validate output |
| `standardize_mol` raises or input parsing returns `None` | Invalid or unsanitizable input | Capture the exception/input index and inspect sanitization deliberately; do not silently accept a partially sanitized structure |
| Stripped wrong fragment | LargestFragmentChooser ambiguity | Manually inspect; consider custom logic |
| Tautomer differs between datasets | Different tautomer rules or toolkit versions | Pin and record the same `TautomerEnumerator` settings and version |
| Unexpected charge distribution with permanent ions | `Uncharger` total-charge policy does not match the intended protonation workflow | Review non-removable and neutralizable sites; choose `force` deliberately and prepare physiological states separately |
| Same InChIKey for apparently different records | Standard-InChI normalization or a rare hash collision | Compare full InChI and standardized structures; InChIKey has no longer form |
| Pipeline slow on large library | Per-molecule Python overhead | Process independent molecules in validated chunks or worker processes; `chembl_structure_pipeline` itself is a per-molecule API |
