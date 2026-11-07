# Drug vs vehicle (drugZ)

```bash
python drugz.py -i experiment.count.txt -o drugz_output.txt -c Veh_r1,Veh_r2 -x Drug_r1,Drug_r2 -p 5
```

`-c` is the vehicle arm, never Day-0. Add `-unpaired` when control and drug replicates do not pair one to one (for example one control against three treatments).

- Output columns: GENE, sumZ, numObs, normZ, then pval/rank/fdr for `synth` (depleted, sensitizers) and `supp` (enriched, suppressors).
- The default `half_window_size=500` needs thousands of guides; a nine-guide demo table fails.
- A second method for consensus: `routes/mle.md` or `routes/rra.md` on the same arms.

Done when `drugz_output.txt` exists.
