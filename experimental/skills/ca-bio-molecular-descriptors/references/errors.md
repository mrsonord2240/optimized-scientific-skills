## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Fingerprint changes between runs | Random seed not set for canonicalization | RDKit Morgan is deterministic; check if input differs (stereo, charges) |
| MACCS bit count != 166 | RDKit MACCS returns 167 bits (bit 0 unused) | Slice `[1:]` if comparing to literature 166-bit |
| Crippen LogP differs from XLogP | Different model | Use `Descriptors.MolLogP` for Crippen; `XLogP3` requires external lib |
| 3D descriptor differs between calls | Different conformer | Set `confId=0` explicitly; or average over ensemble |
| QED returns nan | Charged species or non-standard atom | Standardize (uncharge) before QED |
| Gasteiger charges do not sum to formal charge | Only heavy-atom `_GasteigerCharge` values were summed while hydrogens were implicit | Add explicit H and sum all atom charges, or include each heavy atom's `_GasteigerHCharge` |
| A metal has a neat-looking Gasteiger value | RDKit can emit a finite formal-charge-like value outside the intended organic use case | Reject non-organic elements and use a charge model validated for the element and downstream force field |
| Count-vector similarity differs from bit-vector similarity | Count multiplicities change the generalized Tanimoto calculation | RDKit supports Tanimoto on sparse count vectors; record the vector type and do not compare its threshold directly with a folded-bit threshold |
| MolWt off by ~1 from PubChem | Implicit H counted differently | Use `Descriptors.ExactMolWt` for monoisotopic; PubChem reports average |
