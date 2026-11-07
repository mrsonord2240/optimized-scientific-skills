# Essential genes with BAGEL2

```bash
BAGEL.py fc -i experiment.count.txt -o experiment -c Day0 --min-reads 30
BAGEL.py bf -i experiment.foldchange -o bayes_factor.txt -e CEGv2.txt -n NEGv1.txt \
    -c Day14_r1,Day14_r2,Day14_r3 --seed 42
BAGEL.py pr -i bayes_factor.txt -o pr_curve.txt -e CEGv2.txt -n NEGv1.txt
```

`-c` in `fc` is the baseline column; `-c` in `bf` lists the endpoint replicates. `fc` writes `experiment.foldchange`. `bf` writes `bayes_factor.txt` with columns GENE and BF. `pr` gives the precision-recall used for the CEGv2 PR-AUC gate.

- Always pass a fixed seed, written `--seed`: the short `-s` is overloaded in `bf` (small-sample flag). Unseeded runs are not reproducible and flip hit calls between identical reruns.
- BF above 6 is a hit (FDR near 3%); BF above 3 is FDR near 5%.
- Bootstrapping: add `-b -NB 1000`.
- CEGv2.txt and NEGv1.txt ship in the BAGEL2 clone.

Done when `bayes_factor.txt` exists. Next: `routes/consensus.md` if combining methods.
