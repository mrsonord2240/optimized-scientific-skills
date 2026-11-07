# Acrylamide substitution feature

```bash
python scripts/acrylamide_alpha_substitution.py 'C=CC(=O)N1CCCCC1' > substitution.tsv
```

Replace the example with the input SMILES; the script accepts one or more positional SMILES. None means an invalid input or no match. The count is a structural feature, not a LUMO or reactivity estimate; see `references/design.md`.

Done when requested features and unmatched structures are recorded.
