# Macrocycle embedding

```bash
python scripts/compare_macrocycle_embedding.py --smiles 'O=C1CCCCCCCCCCCNC1' --n-conf 40 --seed 42 > comparison.json
```

Replace the example with the intended macrocycle. The diagnostic compares ETKDGv3 defaults with a torsion opt-out; it is molecule-specific, not proof of ensemble completeness. For the macrocycle helper and limits see `references/macrocycles.md`.

Done when the chosen settings and coverage assessment are recorded.
