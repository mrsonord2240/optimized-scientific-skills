# Score affinity-enrichment data
~~~bash
python scripts/apms_score.py ip.csv --bait Bait1,Bait2,Bait3 --ctrl Ctrl1,Ctrl2,Ctrl3 --out scores.csv
~~~
Supply raw prey-by-replicate intensity or spectral count CSV and negative control IPs. The source defaults require every bait replicate and log2 enrichment >=2 unless absent in all controls; --min-bait-reps changes that tradeoff. Empty controls are reported and all-empty controls fail. This heuristic is not a posterior interaction probability. See [AP-MS notes](../references/affinity_enrichment.md).
Done when controls, detection, sensitivity and ranking are reviewed.
