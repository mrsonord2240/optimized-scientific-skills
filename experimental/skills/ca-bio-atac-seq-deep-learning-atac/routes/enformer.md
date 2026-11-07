# Long-range variant effects

```bash
python scripts/enformer_variant_effect.py genome.fa variants.tsv 12 enformer_scores.tsv
```

Input variants are five headerless fields: chromosome, 1-based position, REF, ALT, ID. Supply actual 0-based output-track indices and document their cell-type match or proxy. This wrapper uses a fixed 196608-base input, forward strand and pseudocount 1; it skips incompatible windows. See `references/method-reference.md`.

Done when selected-track effects, skipped variants and model limitations are reported.
