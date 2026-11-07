# VCF to SPDI
~~~sh
python scripts/query.py vcf 17 43106487 A C result.json
~~~
Supply a valid GRCh38 VCF allele with a 1-based coordinate; the helper builds a 0-based SPDI and asks the source-described canonical-representative endpoint. Verify reference accession and returned warnings, strand and allele. The source helper returns null on unsupported chromosomes and some HTTP errors, so null is not a definitive absence. Source HGVS contextual and SPDI-to-rsID helpers are also available in scripts/dbsnp_lookup.py, but return only the first candidate; inspect the full service response for ambiguity and do not equate contextual with universally canonical representation. See [historical methods](../references/source-methods.md).
**Done when:** the representation, assembly and all unresolved mapping ambiguities are documented.
