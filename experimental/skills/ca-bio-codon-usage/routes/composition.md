# Codon composition

```bash
python scripts/basic_analysis.py
```

This runs a hard-coded illustrative CDS, not a file-input interface. To analyze a supplied CDS, adapt the example and call `count_codons` and `codon_frequencies` in `scripts/codon_utils.py` using the correct genetic-code table and strict validation. Establish orientation and codon position 1 independently; start/stop checks cannot prove either. The example reports counts, frequencies, amino-acid usage and GC123 percentages. Permissive validation discards ambiguous whole codons and trailing partial bases and records them; use it only for explicit cleanup.

**Done when:** the counts, table, frame, sequence length, validation policy and any discarded input are recorded.
