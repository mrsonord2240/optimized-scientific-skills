# Expand site multiplicity
~~~bash
python scripts/expand_sites.py "Phospho (STY)Sites.txt" site_long.tsv
~~~
Source MaxQuant table requires localization, residue, position, protein/gene and bookkeeping columns. This filters class I >=0.75, excludes flagged rows and melts per-run Intensity columns with three underscores and multiplicity 1-3; aggregated columns are excluded. A localization cutoff is not empirical FLR. Review missingness, group ambiguity and per-run columns.
Done when expanded sites, filtering and localization limits are reviewed.
