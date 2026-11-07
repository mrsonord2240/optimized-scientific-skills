# Spectral library search

```bash
python scripts/match_library.py references.mgf queries.mgf --tolerance 0.005 > matches.txt
```

Supply reference and query MGF spectra with identifiable names and precursor m/z. The source applies metadata filters and normalized intensities, ModifiedCosine, score >= 0.7, at least six matched peaks and a 0.02 top-score tie margin. These are screening defaults, not universal identity guarantees. Its printed level is a provisional rule based only on this search, not a comprehensive evidence assessment; inspect candidates, isomers and other evidence. Query names are used as dictionary keys and duplicate or missing names can collide; verify unique identifiers. Spectra with unavailable precursor m/z may fail rather than being excluded. Output redirection overwrites the target.

Done when scores, matched peaks, alternatives, library version and evidence-supported confidence are reviewed.
