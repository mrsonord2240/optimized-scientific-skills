## MAJIQ V3 Differential Analysis

**Goal:** Detect differential LSVs with full posterior distributions over ΔPSI; ideal for complex multi-junction events and heterogeneous cohorts.

**Status: not run.** MAJIQ V3 and VOILA are academic-licence downloads (majiq.biociphers.org) and are not installed here. What follows is taken from the MAJIQ 3.0.11.dev7 documentation (biociphers.bitbucket.io/majiq-docs); check every flag with `majiq <command> --help` before use.

- Pipeline (V3): `majiq build <gff3> <experiments_tsv> <output_dir>` writes `splicegraph.zarr` and one `<experiment>.sj` per sample (options include `--min-experiments`, `--mindenovo`, `--simplify`, `--strandness`); `majiq psi-coverage <splicegraph> <psi_coverage> <sj ...>` prepares coverage (`--minreads`, `--minbins`); `majiq deltapsi` (two groups of replicates, `-grp1/-grp2`, `-n NAME1 NAME2`, `--splicegraph`) and `majiq heterogen` (two groups of independent experiments) quantify.
- V3 replaced V2's `.majiq` files, `-c settings.ini` and SQLite splicegraph; V2-era flags (`--minpos`, `--mem-profile`, `-j`) and `voila view` on V3 output are unverified. The documentation says VOILA "currently only supports MAJIQ v2".
- Use HET for n>=10 vs n>=10 cohort designs (clinical, GTEx-style) and deltapsi for tightly controlled n=3-5. The documentation's statistical guidance is TNOM for n<5 and Wilcoxon for n>5; it does not support a general claim that HET is conservative at n=5-10.
