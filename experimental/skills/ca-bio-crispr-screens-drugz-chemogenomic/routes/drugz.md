# Analyze drug versus vehicle
~~~bash
python drugz.py -i counts.txt -o drugz_output.txt -c Veh_r1,Veh_r2,Veh_r3 -x Drug_r1,Drug_r2,Drug_r3 -p 5
~~~
Run with separately available drugz.py and a TSV whose first columns are guide ID and GENE, followed by sample counts. Use matched vehicle at the same timepoint, not Day 0. Default pseudocount is 5 and half-window 500; small libraries need a smaller window, approximately one quarter of total guides. Inspect sample headers, concordance and guide coverage. Output normZ, fdr_synth and fdr_supp distinguish directions; 0.05 is the source threshold. The source reports a 2019-era API: check local help.
Done when both directional outputs, QC and limitations are reviewed. See [methods](../references/methods.md) and [failure modes](../references/failure-modes.md).
