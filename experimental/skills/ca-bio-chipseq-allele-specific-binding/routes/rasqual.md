# RASQUAL cohort

```bash
python scripts/rasqual_cohort.py prepare --rasqual-source /opt/rasqual --r /opt/R/bin/R --y Y.txt --k K.txt --x X.txt --out cohort-binaries
```

This creates native-double matrices via the source converter in a new directory. Supply rectangular, matching feature-by-sample Y/K tables and sample-by-covariate X table. Prepare a complete feature manifest and indexed phased VCF, then run:

```bash
python scripts/rasqual_cohort.py run --manifest features.tsv --feature-count 2 --covariates 4 --y cohort-binaries/Y.bin --k cohort-binaries/K.bin --x cohort-binaries/X.bin --vcf genotypes.vcf.gz --tabix /usr/bin/tabix --rasqual /opt/rasqual/src/rasqual --output rasqual-cohort.tsv
```

Replace dimensions and paths with actual inputs. The manifest schema, tested source binding, and input conventions are in [methods](../references/methods-and-commands.md). The runner requires one row per feature, checks convergence and output schema, and adds BH adjustment over emitted rows; define that test family before running. Output must not already exist. RASQUAL models mapping bias but sparse data, copy-number imbalance and other confounders remain relevant.

**Done when:** the cohort table, test family, dimensions, phased-genotype provenance and model diagnostics are retained.
