# Motif and signature analyses
~~~bash
python scripts/motif_enrichment.py --sites "Phospho (STY)Sites.txt" --fasta search.fasta --out motif.csv
~~~
This source uses class-I sites, identified-protein background, Fisher and BH; --foreground-ids selects a reviewed subset. For signature input:
~~~bash
python scripts/ptmsea.py write-gct --sites sites.tsv --out sites.gct
~~~
Input seq_window and signed stat must come from an appropriate checked paired-global result for regulatory claims. Duplicate windows are averaged. KSEA, checked ssGSEA execution and result reading are in [kinase](../references/kinase-activity-ksea.md) and [PTM-SEA](../references/ptm-sea.md); never treat proxy results as regulation evidence.
Done when background, contrast direction, prior coverage and outputs are reviewed.

KSEA worker: `scripts/ksea_scores.R`.
