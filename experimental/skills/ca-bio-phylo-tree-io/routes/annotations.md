# Extract tree annotations
~~~bash
python scripts/extract_annotations.py mcc.tree annotations.json --format nexus
~~~
Uses DendroPy metadata extraction and records each node's descendant tip set and annotations per tree. Inspect key names, values, duplicate clade keys and any values stringified for JSON; compare with the original and retain it. This is a side table, not a claim that all annotation semantics round-trip.
Done when needed posterior, HPD and rate fields are verified against the source.
