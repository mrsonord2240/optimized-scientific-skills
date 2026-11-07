# Two factors

The effect of one factor within a level of the other (drug vs vehicle in knockouts): add a combined column to `samples.csv`, such as `group` = `KO_drug`, `KO_vehicle`, `WT_drug`, `WT_vehicle`, then compare two of its levels:

```bash
Rscript scripts/deseq2_de.R counts.csv samples.csv results.csv condition=group ref=KO_vehicle versus=KO_drug
```

Whether the effect differs between levels (is the drug effect different in knockouts): this is the interaction term and needs the code in `references/interaction-designs.md`.

- With `~ genotype + treatment + genotype:treatment`, the coefficient named for treatment is the treatment effect in the reference genotype only. It is not an average over genotypes.
- The combined-column form above avoids that trap. Prefer it whenever the question is about specific pairs.

Done when each requested comparison has its results file.
