# Population frequency

No script; `routes/annotate-vcf.md` attaches the field. Use the grpmax filtering allele frequency (FAF), not a global cutoff.

- FAF is the lower bound of the 95% CI of the grpmax AF (the highest AF among genetic-ancestry groups, formerly popmax). v4 fields: `grpmax`, `AF_grpmax`, `fafmax_faf95_max`, and `fafmax_faf95_max_joint` in the joint exome+genome VCF. Global AF dilutes an allele common in one ancestry.
- A single global cutoff ("AF above 1% is benign") is wrong both ways. The maximum credible AF of a pathogenic allele is per disease (prevalence, genetic and allelic heterogeneity, inheritance, penetrance). Apply BA1/BS1 when FAF exceeds that disease-specific value.
- gnomAD versions: v2.1.1 is 141,456 individuals on GRCh37; v3 is genomes only on GRCh38; v4 is about 730k exomes plus 76k genomes on GRCh38. Comparing AF across versions needs lifting over the variant, not just the coordinate, which can mis-map indels and segmental duplications.
- "Absent" can mean not callable. Check site coverage and the PASS or `AS_FilterStatus` flags. Absence is weaker evidence for an under-represented ancestry.
- "In gnomAD therefore benign" fails for recessive disease (carriers, for example CFTR p.Phe508del), late-onset or reduced-penetrance disease (adult-onset cancer, cardiomyopathy, neurodegeneration), and somatic or clonal-haematopoiesis contamination (low-AF DNMT3A and TET2).
- BA1/BS1/PM2 combining is `variant-calling/clinical-interpretation`.

Done when the frequency field, gnomAD version and build, and the disease-specific threshold used are stated.
