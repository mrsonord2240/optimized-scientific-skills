# Exclude reference genes
~~~bash
python drugz.py -i counts.txt -o drugz_clean.txt -c Veh_r1,Veh_r2 -x Drug_r1,Drug_r2 -r RPS3,RPL11,EIF3A
~~~
The -r argument is a comma-separated list of gene names, not a filename; it removes matching genes before scoring and from output. Compare with an unfiltered run to verify actual matches. See [reference removal](../references/reference-gene-removal.md).
Done when exclusions and their impact are recorded and validated.
