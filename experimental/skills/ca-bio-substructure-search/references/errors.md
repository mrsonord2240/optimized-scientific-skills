## Per-Tool Failure Modes

### PAINS -- false positive on natural product

**Trigger:** Library contains natural products, polyphenols, flavonoids, quinones.

**Mechanism:** RDKit PAINS patterns target rhodanines, catechols, and polyhydroxylated polyphenols -- legitimate scaffolds in natural-product chemistry. RDKit's PAINS implementation does not reproduce every alert class discussed in the original paper; for example, it does not contain a curcumin-specific pattern.

**Symptom:** Library hits flagged as PAINS but trace back to validated natural products with confirmed activity.

**Fix:** Use PAINS as a *flag* not a *delete*. Cross-check flagged compounds for orthogonal-assay confirmation (label-free e.g. SPR, ITC).

### Aromaticity dialect mismatch

**Trigger:** SMARTS pattern with `c` (aromatic) for a heteroatom-rich ring; molecule parsed with different aromaticity model.

**Mechanism:** RDKit, OpenEye, ChemAxon differ on whether furan, thiazole, tropone, etc. are aromatic.

**Symptom:** Same pattern matches in one toolkit, not in another.

**Fix:** Re-canonicalize molecules within RDKit before applying SMARTS. Or use `[#6]:[#6]` instead of `c:c` (explicit element + bond type).

### Tautomer-sensitive pattern miss

**Trigger:** SMARTS targets keto form `C(=O)` but molecule is enol `C(O)=C`.

**Mechanism:** Default canonical form differs by toolkit + standardization choice.

**Symptom:** Known matching molecule reports no match.

**Fix:** Use tautomer-aware match: enumerate tautomers and OR-match. Or canonicalize first via `chemoinformatics/molecular-standardization`. Or expand pattern with `[$(C(=O)),$(C(O)=C)]`.

### Stereochemistry ignored

**Trigger:** SMARTS without `/\@` stereo markers applied to mol with explicit stereo.

**Mechanism:** SMARTS matching is stereo-agnostic by default.

**Symptom:** Wrong stereoisomer is matched as well as right one.

**Fix:** `mol.GetSubstructMatches(pattern, useChirality=True)` to require chirality match.

### Ring closure / fused-ring specificity

**Trigger:** A query must distinguish an isolated benzene ring from a six-membered aromatic ring embedded in a fused system.

**Mechanism:** `c1ccccc1` matches six-membered aromatic cycles and therefore does match benzene cycles within naphthalene. Extra ring-membership or fusion constraints are required to exclude fused systems.

**Symptom:** A nominal "benzene" query returns fused polyaromatics that the project intended to exclude.

**Fix:** Keep `c1ccccc1` when any aromatic six-cycle is desired. When an isolated ring is required, add explicit ring-degree/fusion constraints and test the query against benzene, naphthalene, indole, and representative substituted controls.

### Recursive SMARTS performance

**Trigger:** Deeply nested recursive SMARTS over a large library.

**Mechanism:** Each `[$()]` re-evaluates the inner pattern for every candidate atom.

**Symptom:** Search is materially slower than a flatter query; the magnitude depends on the pattern and library. On a 2,000-molecule reference run, one-level recursion was 1.6x and a deeply nested pattern 3.7x slower.

**Fix:** Flatten recursion where possible; pre-filter with simpler pattern, then re-test with the recursive one.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `Chem.MolFromSmarts` returns None | Invalid SMARTS grammar | Validate with `Chem.MolFromSmarts(smi, mergeHs=False)`; check parens, brackets |
| `[OH]` gives unexpected hydroxyl matches | Query does not state the intended valence/connectivity model | Use `[OX2H]` for neutral alcohol/phenol oxygen or a more specific context-aware pattern |
| Pattern matches but library is "empty" | Mol failed sanitize | Try `Chem.SDMolSupplier(sanitize=False)` then catch errors |
| Multiple matches per molecule | Single-match query expected | `GetSubstructMatch` returns first; `GetSubstructMatches` returns all |
| Match indices but no fragment | Match returns atom indices in pattern order | Map to original mol via `mol.GetAtomWithIdx(i)` |
| PAINS catalog initialization slow | Loading 1000+ patterns on every call | Build catalog once, reuse for batch |
| Stereo SMARTS not matching | `useChirality=False` (default) | `mol.GetSubstructMatches(p, useChirality=True)` |
