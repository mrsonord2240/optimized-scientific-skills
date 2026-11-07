# Published BE Variant-Function Screens

Read when you need the design and results of Hanna 2021 (BRCA1/2, ClinVar-scale) or Cuella-Martin 2021 (86 DDR genes) as methodology references.

## Hanna 2021 BRCA1/2 Variant-Function Screen Methodology

**Hanna et al 2021 *Cell* 184:1064** benchmarked CBE variant scanning at scale, screening 68,526 sgRNAs covering 52,034 ClinVar variants across 3,584 genes, with BRCA1 and BRCA2 as the positive/negative-selection benchmark:

1. Design the CBE library from predicted variant impact (ClinVar annotation), covering each variant with the sgRNAs that install it
2. Run drug-modifier screens (PARPi sensitivity) with vehicle vs drug
3. Score per variant by aggregating over all sgRNAs that install it; cross-check against bystander-controlled sgRNAs

**Standard surrounding practice:** verify editing efficiency at a control timepoint via amplicon sequencing, drop low-efficiency sgRNAs (see the editing-efficiency convention in `screen-analysis.md`), and call sensitizers with a bidirectional method such as drugZ.

**Quantified result:** Recovered known loss-of-function variants in BRCA1 and BRCA2 with high precision, and identified PARP1 variants conferring resistance to PARP inhibitors.

## Cuella-Martin 2021 DDR-Gene Variant Screening

**Cuella-Martin et al 2021 *Cell* 184:1081-1097** screened ~86 DNA-damage-response (DDR) genes (including BRCA1/2) with CBE saturation mutagenesis:

- Saturation CBE design across 86 DDR genes (not BRCA1/2 alone)
- Identified pathogenic/likely-pathogenic variants in critical protein domains
- Combined with biochemical and genetic validation (for example the 53BP1-USP28 interaction surface)
- Demonstrated saturation mutagenesis is feasible at protein-domain scale

**Relationship to Hanna 2021:** the two studies appeared back-to-back in the same *Cell* issue and apply the same CBE variant-scanning strategy to complementary targets -- Hanna benchmarks against ClinVar-annotated variants genome-wide, Cuella-Martin saturates 86 DDR genes. Treat them as complementary methodology references, not as cross-validations of each other.
