# Public example batch annotation
~~~sh
python scripts/query.py batch public_ids.txt annotations.json versions.json
~~~
Use one public identifier per line. This contacts myvariant.info and its metadata service; do not send participant-derived variants without appropriate authorization. The source helper scripts/query_myvariant.py requires myvariant, pandas and requests. Inspect input assembly, allele identity, multi-allelic results and notfound rows; source versions and field paths are historical and unverified. The source describes the service index as GRCh37; check its current contract before preparing HGVS input.
**Done when:** annotations, metadata sidecar, input identifiers and missing or ambiguous results are retained with access date.
