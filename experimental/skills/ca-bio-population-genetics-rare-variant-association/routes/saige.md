# SAIGE-GENE+

```bash
step2_SPAtests.R --bgenFile geno_wes.bgen --bgenFileIndex geno_wes.bgen.bgi --sampleFile samples.txt --groupFile groups.txt --chrom 1 --GMMATmodelFile null.rda --varianceRatioFile null.varianceRatio.txt --annotation_in_groupTest 'lof,missense;lof,missense;lof;synonymous' --maxMAF_in_groupTest 0.0001,0.001,0.01 --is_output_moreDetails TRUE --SAIGEOutputFile gene_tests_chr1.txt
```

This source step-2 example requires an independently fitted compatible null and variance ratio, indexed BGEN, sample file and group file; a complete step-1 pipeline is not supplied. Group variant IDs and annotations must agree. The source describes comma-separated tests and semicolon-joined annotations, and a chromosome-specific LOCO null for each chromosome when LOCO is enabled. Confirm installed help and null compatibility before adapting; no current interface has been verified here. Record the annotation groups, MAF cutoffs, sample order and output per chromosome.

**Done when:** the specified output and model provenance are retained for every intended chromosome with the full test family and masks documented.
