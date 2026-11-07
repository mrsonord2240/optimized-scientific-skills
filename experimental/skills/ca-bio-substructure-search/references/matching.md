## Library Filtering with Multiple Patterns

**Goal:** Reduce a molecule library to those that match all required SMARTS patterns and none of the excluded ones.

**Approach:** Start from the full molecule list, iteratively intersect with each `include` SMARTS using `HasSubstructMatch`, then subtract any molecule matching an `exclude` SMARTS.

```python
def filter_library(mols, include=None, exclude=None):
    keep = list(mols)
    if include:
        for s in include:
            p = compile_smarts(s)
            keep = [m for m in keep if m and m.HasSubstructMatch(p)]
    if exclude:
        for s in exclude:
            p = compile_smarts(s)
            keep = [m for m in keep if m and not m.HasSubstructMatch(p)]
    return keep
```

## Atom Map Indices in SMARTS

Atom maps `[C:1]` track atoms through transformations. Used in reactions (`reaction-enumeration` skill) but also for substructure-based extraction:

```python
# Find amide N with attached aryl
pat = Chem.MolFromSmarts('[CX3:1](=O)[NX3:2][c:3]')
match = mol.GetSubstructMatch(pat)

amide_C, amide_N, aryl_C = match
```
