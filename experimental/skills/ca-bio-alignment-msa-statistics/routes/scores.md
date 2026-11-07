# Substitutions and quality

```bash
python scripts/statistics.py scores alignment.fasta scores.json
```

The result contains raw unordered mismatch counts, flat-gap alignment score and, for nucleotide data, transitions and transversions. These are not a log-odds substitution matrix. Use `references/substitution-counts.md` for pairwise substitutions, `references/alignment-quality-scores.md` for BLOSUM62 sum-of-pairs, and `references/information-content-pssm.md` for PSSMs. Alphabet-bound scores omit unsupported letters; report omissions.

For the reference code imports and normalized alignment setup use `references/setup.md`.

Done when the requested score is computed and its definition reported.
