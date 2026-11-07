# MOL2 and PDBQT conversion

```bash
obabel input.mol2 -O output.sdf
```

Open Babel is a route for MOL2 dialects and PDBQT. Inspect stereo and aromaticity after conversion; use RDKit for registry InChI unless the Open Babel build has confirmed support. See `references/conversion.md` and `references/failures.md`.

Done when converted structures and any failed or altered records are reviewed.
