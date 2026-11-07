# Long sequences or large searches

```bash
python scripts/align.py sequences.fasta alignment.json --mode local --type dna
```

For large pairs, edit distance, high-throughput or read mapping select a library in `references/library-selection.md`; check SIMD saturation. For query-versus-database or very low identity use `references/when-not-appropriate.md` and canonical `bio-alignment-structural` when structures are available. For three or more sequences use `bio-alignment-multiple`; codon-aware alignment belongs there as well.

Done when the requested comparison or search uses a method suited to its input and scoring requirements.
