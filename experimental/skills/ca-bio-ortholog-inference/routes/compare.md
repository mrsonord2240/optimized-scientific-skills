# Compare resources

```bash
python examples/cross_resource.py
```

This source demonstration compares TP53 across Compara, OMA and a verified OrthoDB group; adapt the input IDs and group search before applying it to another gene. It reports resource failures separately from empty results, but does not normalize identifiers or calculate an intersection. Normalize identifiers to a shared namespace, compare calls and inspect disagreements manually. Confidence semantics differ by resource and cannot be averaged; see [interpretation](../references/interpretation.md). For custom species outside databases use a de novo method; this source does not implement one.

Done when resource versions, identifier mapping, overlap, disagreement and unavailable resources are documented.
