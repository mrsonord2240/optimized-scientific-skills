# Variant consequences and LD

```bash
python scripts/query.py /vep/human/id/rs55794205 variant.json
```

Choose actual species, variant ID, HGVS or region and allele; match coordinates to the assembly. VEP and pairwise LD endpoints are listed in `references/endpoints.md`. For large variant sets use a local VEP workflow rather than long per-variant REST loops. Returned predictions are annotations and require interpretation in context.

Done when the requested response and assembly, release and input variant are recorded.
