# Evaluate a pairwise score

```bash
python scripts/align.py sequences.fasta alignment.json --mode local --type protein
```

Do not compare raw scores across scoring systems or infer significance from an alignment alone. For calibrated bit scores, E-values, composition masking and source empirical-shuffle implementation read `references/significance.md`. For DNA use the dinucleotide-preserving shuffle where appropriate. Sequence-only evidence becomes unreliable at low identity; `references/when-not-appropriate.md` lists profile and structural alternatives.

Done when the score is interpreted with its scoring scheme and appropriate significance or limitations.
