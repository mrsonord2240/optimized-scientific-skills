# Sequence Filtering

Read when the request is to subset sequences by ID pattern, gap content or uniqueness. Script: `scripts/filter_sequences.py`.

## Sequence Filtering

**Goal:** Subset an alignment to retain only sequences matching specific criteria (ID pattern, gap content, uniqueness).

**Approach:** Iterate over alignment records, apply filter conditions, and reconstruct a new MultipleSeqAlignment from matching records. A filter that would remove every sequence raises `ValueError` instead of returning an empty alignment.

```bash
python scripts/filter_sequences.py in.fasta out.fasta --id-pattern '^species_' --max-gap-fraction 0.1 --dedup
```

```python
import sys; sys.path.insert(0, 'scripts')  # run from this Skill's directory
from filter_sequences import filter_by_id, filter_by_gap_content, remove_duplicates
```

`remove_duplicates` compares normalised rows (`AC-GT`, `AC.GT` and `ac-gt` are one sequence) and keeps the original records; the kept rows keep the input's annotations.
