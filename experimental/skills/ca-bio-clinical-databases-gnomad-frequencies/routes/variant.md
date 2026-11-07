# Public variant frequency lookup
~~~sh
python scripts/query.py variant 17 43106487 A C result.json
~~~
This sends a public GRCh38 example to the source-described gnomAD v4 GraphQL endpoint. The helper scripts/gnomad_query.py imports requests, myvariant and pandas. Verify installed APIs and actual schema; this draft makes no live calls. A not-found response can also reflect a wrong build or allele and does not prove absence. Review filters, allele number, release, ancestry and FAF95 missingness. Never send patient-derived variants to a public API without appropriate authorization.
**Done when:** raw response and extracted frequency provenance are retained, with uncertainty documented.
