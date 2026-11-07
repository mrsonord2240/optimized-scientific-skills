# Public gene constraint lookup
~~~sh
python scripts/query.py gene SCN2A result.json
~~~
The source query obtains observed/expected metrics for a gene where available. Constraint and transcript depend on release; compare within an appropriate release and inspect the actual API response. The source includes contradictory v4 first-decile thresholds, so this draft does not calculate a decile label automatically. Do not infer patient-level pathogenicity from gene constraint.
**Done when:** metrics, transcript, release and missingness are documented for qualified review.
