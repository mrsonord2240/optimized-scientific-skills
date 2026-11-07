# Bridge TMT plexes
~~~bash
python scripts/normalize_matrix.py irs --plex plex_a.csv --ref A_ref --plex plex_b.csv --ref B_ref --out-dir bridged
~~~
Each CSV is protein by channel with a pooled reference; sample-loading normalization precedes IRS. Proteins with a zero or missing reference in any plex become unbridged NaNs; review them and check non-reference channels. Reporter extraction and lot-specific isotope impurity correction remain source R snippets in [TMT notes](../references/tmt_isobaric.md); this helper does not correct impurities or ratio compression.
Done when reference alignment, unbridged proteins and within-plex normalization are reviewed.
