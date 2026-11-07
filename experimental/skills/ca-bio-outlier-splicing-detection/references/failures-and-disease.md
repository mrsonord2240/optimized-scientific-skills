# Per-Tool Failure Modes


### FRASER 2: Tissue Mismatch

**Trigger:** Patient sample from a different tissue than the controls. **Mechanism:** the fit learns tissue-specific patterns; a mismatched patient becomes a global outlier or is absorbed. **Symptom:** many implausible calls (AE) or none (PCA); see Tissue Choice for the check. **Fix:** use only controls from the patient's tissue.

### OUTRIDER: Few Controls

**Trigger:** Cohort <50 samples. **Symptom:** no convergence warning; `Warning: No significant events` (measured n=8-30, 0 false calls). **Fix:** recruit >=50 matched controls; a z-score on log-CPM against the controls is a cruder fallback. Do not read an empty result as "no expression outliers".

### LeafcutterMD: Cluster Count Limits

**Trigger:** Very few clusters in the patient sample (low coverage or filtered out). **Mechanism:** the Dirichlet-multinomial fit is unstable with few observations. **Symptom:** inflated or deflated p-values; few significant calls. **Fix:** increase coverage; relax filtering (`-m 10` instead of 50); or switch to FRASER2.

### DROP: Snakemake Pipeline Failures

**Trigger:** Missing dependencies or incompatible R/Bioconductor versions; a killed run leaves a lock. **Symptom:** a step fails partway with a cryptic R error. **Fix:** install DROP in its own bioconda environment (see DROP Pipeline); after a killed run `snakemake --unlock`, then resume.

## Reconciliation: When Outlier Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| FRASER2 sig, OUTRIDER not | Splicing change without expression collapse | Standard splicing outlier; report |
| OUTRIDER sig, FRASER2 not | Expression LoF without splicing change | Likely promoter / regulatory; not splicing |
| Both sig at same gene | LoF allele triggering NMD on splicing-disrupted transcript | Strong combined evidence; expect downstream |
| LeafcutterMD sig, FRASER2 not | Novel cryptic event not in annotation | High-priority novel finding; investigate |
| All tools null but biology suggests change | Underpowered cohort or wrong tissue | Verify gene expression in tissue; recruit larger cohort |

## Disease Context

Use disease biology to choose the tissue and candidate follow-up, but do not treat a familiar signature as an expected result or a diagnostic conclusion. Published examples include TDP-43-associated cryptic exons in ALS/FTD and SF3B1-associated alternative 3' splice-site use (references below). For every case, first verify that the candidate gene is expressed in the queried tissue (for example with GTEx) and interpret a negative RNA result as inconclusive when it is not.
