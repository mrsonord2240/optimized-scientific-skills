# Condition and classification review
A VCF or VCV summary is not a substitute for condition-specific RCV and individual SCV evidence. The source describes RCV XML parsing, VCEP review, conflict evaluation, somatic clinical impact and oncogenicity, but does not provide a complete condition-specific XML parser. Adapt a version-verified parser and retain distinct interpretation types, submitters, dates and unresolved conflicts. Review the [source methods](../references/source-methods.md) and [usage guidance](../references/usage-guide.md); historical release and schema claims are unverified.
For local annotation, the source gives:
~~~sh
bcftools annotate -a clinvar.vcf.gz -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVC,INFO/CLNHGVS,INFO/CLNSIGCONF input.vcf.gz -Oz -o annotated.vcf.gz
~~~
Use matching assemblies and contigs and inspect the VCF header for fields and current definitions. Related shelf skills: bio-clinical-databases-acmg-classification, bio-clinical-databases-gnomad-frequencies and bio-clinical-databases-variant-prioritization.
**Done when:** classifications, review status, condition, conflict, source version and limits are documented for qualified review.
