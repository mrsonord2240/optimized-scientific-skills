# Cumulative incidence with competing events
~~~sh
python scripts/cif.py outcomes.csv time cause 1 cif.csv
~~~
The cause column uses 0 for censoring and positive integer event codes; 1 is the chosen event in this example. Aalen-Johansen estimates the cumulative incidence under competing risks. Record time origin, follow-up and censoring assumptions; source notes lifelines may jitter tied event times; this wrapper sets seed 0 and requires review of tie handling. This is a population research estimate, not an individual prognosis.
**Done when:** the curve, event definition, competing causes and uncertainty are reviewed.
