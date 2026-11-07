# Import MaxQuant LFQ
~~~bash
python scripts/load_maxquant.py proteinGroups.txt log2_lfq.tsv
~~~
The source reader filters present Reverse, Potential contaminant and Only identified by site flags, retains first semicolon IDs, converts zero to NaN and log2-transforms LFQ. It requires Protein IDs and at least one flag column; cross-check REV__/CON__ IDs and audit absent flag columns. The first ID does not resolve all protein-group ambiguity. Review contaminants before removal if they are a target. Prints row counts and an abundance/missingness correlation; correlation alone cannot establish a missingness mechanism.
For corrected TMT reporter channels the source supplies a conditional code example in [quantitation](../references/quantitation.md), not a bundled CLI. R users may run scripts/load_maxquant_qfeatures.R on the input for a QFeatures summary; its process exit must be zero.
Done when filtering counts, quant family, missingness and output are reviewed.
