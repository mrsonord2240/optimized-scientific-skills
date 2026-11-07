# Phylogenetic alignment

```bash
clipkit input.fasta -m smart-gap --log -o trimmed.fasta
```

For gap- or entropy-based filtering, retention below 0.6 signals an overly aggressive setting: inspect and compare with less filtering or the original. This heuristic does not apply to kpi/kpic modes, which deliberately remove singletons; compare topology and branch lengths, including outgroups, instead. `kpic-smart-gap` is a topology-focused choice for balanced samples, not branch-length, rate or dating work. Short genes may be better left untrimmed. See `references/clipkit.md` for modes and `references/mapping.md` for output coordinates.

Done when alignment, log and retention are checked and tree sensitivity to trimming is reported.
