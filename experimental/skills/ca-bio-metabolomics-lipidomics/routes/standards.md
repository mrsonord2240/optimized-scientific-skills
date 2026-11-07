# Class internal-standard normalization

```bash
Rscript scripts/istd_normalize.R out_dir annotation.csv skyline1.csv
```

Supply raw Skyline export and sample annotation; multiple Skyline files can follow. The script warns and drops unclassified lipids, stops if a recognized class lacks an internal standard, then writes log2 normalized values and a class plot. Check dropped names and standards against spike metadata, extraction timing and input format. Supplying only an output directory invokes bundled demo data, not the user's data. Existing outputs may be overwritten. Ratios are within-class and do not establish cross-class molar comparisons without response-factor calibration.

Done when all retained classes have appropriate standards and exclusions, normalization and interpretation limits are documented.
