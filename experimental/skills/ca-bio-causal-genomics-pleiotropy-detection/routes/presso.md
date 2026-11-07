# MR-PRESSO
~~~sh
Rscript scripts/presso.R harmonized.tsv presso.rds 5000
~~~
The optional final argument sets simulations (default 5000); seed is 42. Supply harmonized beta and standard-error columns and at least four retained instruments. Use enough simulations for the number of instruments and desired significance. A non-significant global test does not rule out correlated horizontal pleiotropy. Review outliers and distortion in context.
**Done when:** the RDS contains the global and requested outlier and distortion results and their limitations are documented.
