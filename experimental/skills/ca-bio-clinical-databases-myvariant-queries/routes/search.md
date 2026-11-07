# Public gene search
~~~sh
python scripts/query.py gene BRCA1 hits.json 500 2
~~~
This contacts myvariant.info. The optional final arguments are result limit (default 500) and minimum review-star filter (default 2). Compare total with returned hits; a truncated result is not exhaustive. Source query terms may match different RCV entries, so review condition and submission-level evidence individually. The source also provides a GRCh37 region search in scripts/query_myvariant.py; adapt and verify field paths and coordinates first.
**Done when:** query, total, returned records, source versions and any truncation are documented.
