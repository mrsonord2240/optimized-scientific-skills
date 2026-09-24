# Interaction designs in DESeq2

Interaction design with the canonical trap:

```r
design(dds) <- ~ genotype + treatment + genotype:treatment
dds <- DESeq(dds)
resultsNames(dds)
# "Intercept" "genotype_KO_vs_WT" "treatment_drug_vs_vehicle" "genotypeKO.treatmentdrug"
```

- `results(dds, name='treatment_drug_vs_vehicle')` returns the drug effect IN THE WT REFERENCE only, NOT a marginal average. This is the single most common misinterpretation.
- `results(dds, name='genotypeKO.treatmentdrug')` returns the DIFFERENCE in drug effect between KO and WT.
- Drug effect in KO requires summing: `results(dds, contrast=list(c('treatment_drug_vs_vehicle','genotypeKO.treatmentdrug')))`.

Cleaner alternative for interactions when many pairwise contrasts are needed: combined factor + `~ 0 + group`.

```r
dds$group <- factor(paste(dds$genotype, dds$treatment, sep = '_'))
design(dds) <- ~ 0 + group
dds <- DESeq(dds)
res_drug_in_ko <- results(dds, contrast = c('group', 'KO_drug', 'KO_vehicle'))
```
