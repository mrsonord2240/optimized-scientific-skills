# Pairwise identity

```bash
python scripts/statistics.py identity alignment.fasta identity.json --method pid4
```

Report the denominator and alignment length. PID1 uses aligned pairs plus internal gaps; PID2 uses aligned residue pairs; PID3 uses the shorter ungapped length; PID4 uses mean ungapped length. PID4 is the source default for full-length evolutionary comparisons; PID2 describes an aligned domain or fragment. Undefined values are null. Exact pairwise calculation scales quadratically in sequence count; consider k-mer estimation for very large datasets.

Done when the requested matrix or pairwise values are reported with the chosen definition.
