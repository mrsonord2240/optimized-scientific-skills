# Merge and inspect transition libraries
~~~bash
python scripts/library_tools.py merge merged.tsv library_a.tsv library_b.tsv
~~~
The source combines inputs and deduplicates by modified sequence, precursor charge, fragment type, series number and fragment charge, preferring rows sorted by precursor-total intensity. It reports precursor and protein counts and mean transitions per precursor. Check input units, intensity comparability and conflicting duplicate metadata before accepting the result; see [QC notes](../references/library-qc-merge.md).
Done when preserved transitions, metadata conflicts and summary are reviewed.
