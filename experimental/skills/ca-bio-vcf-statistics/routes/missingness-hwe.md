# Missingness, call rate and HWE

```bash
bash scripts/missing_hwe.sh input.vcf.gz qc_out [controls.vcf.gz]
```

Writes `sample_miss.imiss`, `site_miss.lmiss` and exact-test `hwe.hwe` (from `controls.vcf.gz` when given, else the input) and prints the 10 worst samples by `F_MISS`.

- Genotype-level filters (`./.` where GQ<20 or DP<8) come first; low-quality genotypes drive both missingness and HWE.
- Sample and variant missingness are coupled: drop the worst samples, re-run, drop the worst variants, repeat. Do not apply both thresholds in one pass. Conventions drop samples below ~95-98% call rate and variants below ~95-99%.
- Filter on excess heterozygosity only. Het deficit is often real (Wahlund effect, inbreeding, selection, null allele).
- Run HWE within one ancestry-homogeneous group; a mixed cohort induces a het deficit that mimics genotyping error.
- Use the exact test (`vcftools --hardy`), not chi-square, which is anticonservative for rare variants.
- In case/control studies compute HWE in controls only; a real strong-effect association deviates in cases.
- GATK `ExcessHet` (phred; larger is more suspect), default cutoff 54.69, tune to cohort size:

```bash
bcftools query -f '%CHROM\t%POS\t%INFO/ExcessHet\n' input.vcf.gz | awk '$3>54.69'
```

Done when missingness and HWE files exist and the filters you applied name their order.
