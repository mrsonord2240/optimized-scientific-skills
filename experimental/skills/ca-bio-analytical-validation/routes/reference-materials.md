# Choosing a reference material

No script. Record the exact material, dilution, fragmentation state, truth set and commutability caveat in the validation plan.

| Validating | Use | Identity |
|------------|-----|----------|
| Somatic variant callers (paired tumour-normal) | HCC1395 tumour (ATCC CRL-2324) with HCC1395BL matched normal (ATCC CRL-2325) and the SEQC2 somatic call set v1.2 | SRA SRP162370 |
| Oncopanel or ctDNA workflows at low allele frequency | SEQC2 Sample A / Sample B-derived dilutions, for example Df and Ef fragmented to 130-170 bp | BioProject PRJNA677999 |

- Sample A is an equal-mass pool of ten UHRR cancer-cell-line DNAs, not HCC1395. Sample B is Agilent male reference DNA 5190-8848.
- Contrived reference standards are not commutable with patient plasma; say so.
- Dilution-series design for the LoD fit: `routes/lod95-probit.md`.

Done when the plan names one material per question and states its caveat.
