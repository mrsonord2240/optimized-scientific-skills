# Gap statistics

Read when summarising the gap distribution of an alignment (per column, per sequence, gap-free columns). Assumes the normalised `alignment` and the imports from "Required Import and Normalisation" in `SKILL.md`.

## Gap Statistics

**Goal:** Summarize gap distribution across the alignment to assess alignment quality and identify problematic regions.

**Approach:** Calculate gap fractions per column and aggregate statistics including total gaps, gap-free columns, and gappiest sequence/column.

### Gap Fraction Per Column
```python
def gap_profile(alignment):
    profile = []
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        gap_fraction = column.count('-') / len(alignment)
        profile.append(gap_fraction)
    return profile

gaps = gap_profile(alignment)
avg_gaps = sum(gaps) / len(gaps)
print(f'Average gap fraction: {avg_gaps*100:.1f}%')
```

### Gap Statistics Summary

```python
def gap_statistics(alignment):
    # total gaps, gap fraction, gappiest sequence and column indices, gap-free column count
    ...
```

Full implementation: `examples/gap_statistics.py` (`gap_statistics(alignment)` returns a dict).
