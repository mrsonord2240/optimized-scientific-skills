# Resolve a public rsID
~~~sh
python scripts/query.py resolve rs630496 result.json
~~~
This contacts the source-described NCBI Variation Services endpoint. The helper follows merge hops with a cycle check and ten-hop limit, returning resolved, not_found, withdrawn, orphan or an error. A 404 does not prove withdrawal. Inspect the full record and all matching alleles before downstream use. Service paths, schemas and rate limits have not been reverified. Do not send patient-derived identifiers to public services without appropriate authorization.
**Done when:** the result, chain, access date and unresolved status are recorded.
