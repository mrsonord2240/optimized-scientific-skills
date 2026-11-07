# Pooled bridge layout
~~~sh
Rscript scripts/bridge_layout.R samples.csv 4 16 16 condition,site layout.csv 10000 1
~~~
Inputs are sample CSV, plex count, channels per plex, reserved channel, comma-separated balance features, output, optional iterations and seed. The helper reserves one channel per plex, checks capacity and calls the same conservative coverage check as [assignment](assignment.md). The pooled reference is not a biological replicate. Randomize injection order independently and verify actual designit API before use.
**Done when:** reserved channels remain unassigned, sample constraints and full-model estimability are reviewed, and the layout is recorded.
