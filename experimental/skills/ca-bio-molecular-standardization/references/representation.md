## Salt Stripping Edge Cases

| Salt form | Action | Example |
|-----------|--------|---------|
| Mono-salt | Strip counter-ion | `[Na+].CC(=O)[O-]` -> `CC(=O)O` |
| Di-salt | Strip both | `[Na+].[Na+].CC(=O)[O-].CC(=O)[O-]` -> `CC(=O)O` |
| Mixed salt | Largest organic fragment | `CCO.CC(=O)O` -> `CCO` (or CC(=O)O depending on rule) |
| Co-crystal | Hardest case | `CC(=O)O.CCOC(C)=O` -- both organic; default returns largest |
| Hydrate | Strip waters | `CC(=O)O.O` -> `CC(=O)O` |
| Solvate | Strip solvents | `CC(=O)O.CO` -> `CC(=O)O` |
| Quaternary ammonium | Preserve charge | `[N+](C)(C)(C)C` (permanent charge; do NOT neutralize) |

**`LargestFragmentChooser(preferOrganic=True)`** prefers organic fragments over inorganic counter-ions even if smaller; for co-crystals, default rule picks largest organic fragment.

## Tautomer Canonicalization (debated)

Tautomer canonicalization is the most controversial standardization step. There is no universally-correct canonical tautomer for many drug-like molecules.

| Tautomer pair | Why the policy matters |
|---------------|------------------------|
| Keto/enol | Canonicalization can select a representation different from the experimentally relevant bound or solution form |
| Lactam/lactim | Heterocycle scoring rules and toolkit versions may choose different representatives |
| Amidine/iminol | Proton placement changes donor/acceptor annotations and downstream matching |
| Phenol/keto (e.g., naphthol/naphthalenone) | Aromaticity and functional-group perception can change with the selected representation |
| 2H-pyrazole / 1H-pyrazole | Nitrogen identity and donor/acceptor assignments depend on proton placement |

Treat the enumerator's canonical result as a reproducible representation chosen by its configured scoring rules, not as a prediction of the dominant tautomer in vivo. Record the RDKit version and any custom transforms or scoring changes.

**Practical rules:**
- Always apply consistent canonicalization across train + test for ML
- For prospective prediction, predict for both tautomers if disagreement could matter
- For library deduplication, canonical tautomer is the standard answer
- For docking, use an ionization-aware preparation workflow. For Open Babel, the documented CLI is `obabel input.sdf -O output.sdf -p 7.4`; validate generated states because its rule-based protonation is not a substitute for project-specific pKa analysis.

```python
from rdkit.Chem.MolStandardize import rdMolStandardize

def canonical_tautomer(smi):
    mol = Chem.MolFromSmiles(smi)
    enumerator = rdMolStandardize.TautomerEnumerator()
    canon = enumerator.Canonicalize(mol)
    return Chem.MolToSmiles(canon)
```

## Stereochemistry Standardization

```python
from rdkit import Chem

def standardize_stereo(mol, remove_undefined=False):
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if remove_undefined:
        Chem.RemoveStereochemistry(mol)
    return mol
```

**Cases:**
- Explicit stereo with `@` / `\` / `/` -> preserved
- Wedge bonds in SDF -> re-perceived from 3D coords if present
- Ambiguous stereo (no markers) -> left as-is, marked as undefined
- Racemic (explicit "rac") -> keep as racemate

For ML, remove stereochemistry only when the endpoint, data curation, and model representation justify treating stereoisomers as equivalent; record that policy and test its effect. For docking and FEP, preserve the intended stereoisomer and reject unintended stereo changes.
