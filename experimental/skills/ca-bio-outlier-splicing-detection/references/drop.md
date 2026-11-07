# DROP Pipeline (Integrated Workflow)


**Goal:** Run FRASER2 + OUTRIDER + monoallelic expression in one Snakemake pipeline.

**Approach:** DROP is distributed via **bioconda** (not PyPI). Install in a dedicated environment, initialise an empty project directory, fill in the sample annotation and `config.yaml`; the pipeline handles counting, autoencoding and reporting.

```bash
# Install via bioconda (DROP is not on PyPI)
mamba create -n drop_env -c conda-forge -c bioconda drop --override-channels
conda activate drop_env

mkdir my_diagnostic_run && cd my_diagnostic_run
drop init          # takes no project-name argument; run it in the empty directory

# Edit config.yaml (checked on DROP 1.6.1): fill sampleAnnotation, geneAnnotation, genome, root, htmlOutputPath;
#  aberrantSplicing: run: true   (implementation: PCA, FRASER_version: "FRASER2", padjCutoff 0.1, deltaPsiCutoff 0.1)
#  aberrantExpression: run: true (implementation: autoencoder)
#  mae: run: true                (needs RNA BAMs matched to a VCF)

snakemake aberrantSplicing --cores 16    # or aberrantExpression / mae; no --use-conda: DROP's rules have no conda: directives
```

`drop demo` downloads a public 10-sample demo project; its `snakemake aberrantSplicing` took about 90 minutes to first `results.tsv` on 10 cores and reported no significant outlier (min padjust 1), as expected for 10 samples. Allow longer for all targets and a cold cache. The **MAE module** tests allelic imbalance on heterozygous SNPs called from RNA-seq with a negative-binomial test (not a z-score); useful for monoallelic LoF that splicing/expression outliers miss.
