# ClinVar research lookup
~~~sh
python scripts/lookup.py local clinvar.vcf.gz 17 43106487 A C result.json
~~~
The source helper scripts/clinvar_query.py imports requests, pandas and cyvcf2 for every mode. Use a matching indexed local ClinVar VCF and the correct genome assembly, normalized alleles and contig names. The source example expects unprefixed contigs and GRCh38. This returns variant-level VCF INFO, not a condition-specific classification; missing results do not imply benignity. For a public VariationID, the source also provides:
~~~sh
python scripts/lookup.py summary 17661 result.json
~~~
The summary command contacts NCBI. Its source-described schema and review-status mapping are unverified here; inspect actual responses and service policies. Never send patient-derived identifiers to public services without appropriate authorization.
**Done when:** output, source release, assembly and query identity are retained and any missing or conflicting evidence is reviewed.
