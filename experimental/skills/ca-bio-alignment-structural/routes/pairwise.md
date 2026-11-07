# Pairwise structural alignment

```bash
TMalign mobile.pdb reference.pdb -o superposed > alignment.txt
```

The prefix produces `superposed.pdb`. Verify data were emitted: short or unparsable input can exit successfully without an alignment. For known, co-numbered atom correspondence use `python scripts/biopython_superimposer.py reference.pdb mobile.pdb`; it writes `mobile_superposed.pdb` and rejects poor correspondences. For complexes use USalign as described in `references/pairwise.md`. Report both TM normalizations and the aligned length; use the smaller TM for whole-chain comparison, and avoid fold calls from TM alone below 60 residues. For GDT-TS use `TMscore model.pdb native.pdb -seq`; residue numbering alone can mispair.

Done when the structure and metrics are written and their correspondence and normalization are stated.
