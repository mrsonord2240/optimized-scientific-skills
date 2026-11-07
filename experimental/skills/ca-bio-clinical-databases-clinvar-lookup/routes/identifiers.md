# Allele Registry resolution
~~~sh
python scripts/lookup.py allele 'NC_000017.11:g.43106487A>C' allele.json
~~~
This contacts the ClinGen Allele Registry with a public HGVS example. Record the reference accession and returned record; mismatch or HTTP errors require review, not a benign interpretation. The example helper also supports batch resolution, but first check consent, governance, rate limits and the current service contract. A CA ID, VariationID, VCV, RCV and SCV are different identifiers; do not assume a mapping without examining the record.
**Done when:** the returned public record and provenance are saved or a missing/error response is documented.
