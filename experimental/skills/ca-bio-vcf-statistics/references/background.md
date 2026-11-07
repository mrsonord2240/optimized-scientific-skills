# Background

Not needed to run any route.

## Why no metric has a universal pass value

Each metric's expected range depends on the assay (WGS vs WES), the sample's ancestry and the cohort. QC is reading a deviation for its mechanistic meaning, not threshold-checking. Order matters: genotype-level filters before cohort missingness, call rate and HWE, because low-quality genotypes in the matrix drive spurious missingness and HWE deviation. Sample identity QC is a graph problem: build the all-pairs relatedness matrix and reconcile it against the manifest before any association or burden analysis, since one undetected swap or contaminated sample can fabricate or erase a genome-wide-significant hit.

## Ti/Tv

Transitions are favored over transversions because CpG deamination (methylated C to T) is the most common vertebrate point mutation and is a transition, and because transitions are more often synonymous. A random error spectrum gives Ti/Tv ~0.5, so false positives pull the ratio down. WES runs higher than WGS because coding regions are CpG- and constraint-enriched and codon degeneracy favors transitions.

## Het/hom and ancestry

The het:non-ref-hom ratio reflects heterozygosity relative to the reference. African-ancestry genomes diverge more from GRCh (a mostly European-ancestry assembly), so a single cutoff tuned on EUR would flag every AFR sample.

## Novel versus known

Real common variants are already catalogued, so novel common variants are almost always artifacts; novel rare or singleton variants are expected.

## HWE

Heterozygote excess is the signature of a mapping artifact: a duplicated or collapsed region piles reads from two paralogous copies onto one locus. GATK `ExcessHet`, `InbreedingCoeff` and gnomAD filter the excess side only. The 54.69 cutoff is GATK's default (phred-scaled p ~3.4e-6, about the 1000-sample z=-4.5 boundary).

## Contamination

Foreign reads add minor-allele support at true hom sites and skew true hets, shifting the het allele-balance distribution off 0.5.
