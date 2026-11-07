# Assess temporal structure
~~~bash
python scripts/root_to_tip.py tree.nwk dates.tsv
~~~
Dates TSV contains tip and decimal date; source searches edge-midpoint rootings for maximum R-squared, which is exploratory and does not account for non-independent tips. Positive slope alone is insufficient. For a source date-randomization approximation:
~~~bash
python scripts/date_randomization.py alignment.fasta tree.nwk dates.tsv 20 42
~~~
Run in a fresh dedicated directory: the source uses fixed prefixes and IQ-TREE -redo, which can overwrite matching outputs. It compares the real point estimate against the shuffled point-estimate range, not confidence-interval overlap; it also parses output without checking subprocess return code. Review failures, dependence, precision and formal test design before claiming temporal signal.
Done when data provenance, rooting, randomization assumptions and uncertainty are reviewed.
