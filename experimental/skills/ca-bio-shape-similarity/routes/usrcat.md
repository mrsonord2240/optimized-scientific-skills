# USRCAT search

```bash
python scripts/usrcat_search.py query.smi library.smi usrcat_scores.csv
```

Each file contains one SMILES per line; the first query is used. The wrapper embeds one conformer per molecule with seed 42 and records failures. USRCAT is a coarse alignment-free prefilter; see `references/usrcat.md` and `references/selection.md`.

Done when input-ordered scores and failed records are saved.
