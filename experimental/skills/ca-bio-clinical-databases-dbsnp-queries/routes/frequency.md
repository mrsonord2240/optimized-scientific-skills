# ALFA frequencies
~~~sh
python scripts/query.py alfa rs6025 frequencies.json
~~~
The helper queries source-described ALFA allele counts by population; check the returned placement and allele. ALFA and gnomAD differ in ascertainment and population composition; raw ALFA AF is not FAF95. The source batch helper collapses results and requires myvariant, pandas and requests; validate per-allele joins and failures rather than treating a missing record as zero. Genome-wide local JSON streaming is described but no complete reusable parser is supplied. Consult [usage guidance](../references/usage-guide.md).
**Done when:** release, population, allele, counts and missingness are retained for research review.
Related shelf skills: bio-clinical-databases-clinvar-lookup, bio-clinical-databases-gnomad-frequencies and bio-clinical-databases-myvariant-queries.
