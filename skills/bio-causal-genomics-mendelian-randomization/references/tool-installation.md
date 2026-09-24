# Tool installation

Read when a package is missing or `install.packages` fails. Verbatim from SKILL.md "Tool Installation Notes".

## Tool Installation Notes

```r
# CRAN-stable
install.packages(c('remotes', 'MendelianRandomization', 'MVMR', 'coloc', 'simex'))

# GitHub-only or recently archived
remotes::install_github('MRCIEU/TwoSampleMR')          # primary orchestrator
remotes::install_github('MRCIEU/ieugwasr')             # OpenGWAS client + local clumping
remotes::install_github('rondolab/MR-PRESSO')          # never on CRAN
remotes::install_github('qingyuanzhao/mr.raps')        # CRAN-archived 2025-03-01
remotes::install_github('jean997/cause')               # depends on mixsqp; suggests Rfast
remotes::install_github('cnfoley/mrclust')             # heterogeneity clusters
remotes::install_github('LizaDarrous/lhcMR')           # bidirectional + heritable confounder
remotes::install_github('HDTian/DRMR')                 # doubly-ranked stratification
remotes::install_github('n-mounier/MRlap')             # joint overlap + winner's-curse + weak-IV correction
```

`TwoSampleMR::mr_raps()` is a thin wrapper that calls `mr.raps::mr.raps()` under the hood; the GitHub `mr.raps` install above is therefore required. The `MendelianRandomization` package does NOT export `mr_raps()` (verify with `ls('package:MendelianRandomization')`); only TwoSampleMR offers a MR-RAPS entry point.

**Local clumping (verified 2026-09-21):** `install.packages('remotes'); remotes::install_github('MRCIEU/genetics.binaRies')` -- a small package (~1s install; its bundled `plink.exe` is fetched on first `get_plink_binary()` call, ~28 MB, cached after). Do NOT install plink2 for this: `ieugwasr::ld_clump_local()` shells out to plink **1.9** syntax (`--clump <file> --clump-p1 --clump-r2 --clump-kb`, reading the `.clumped` output file plink1.9 writes); plink2's `--clump` takes a different report format and writes a different output filename, so it does not work as a drop-in. `genetics.binaRies::get_plink_binary()` is what `scripts/twosample_workflow.R` calls by default and is the one that actually works. For the reference bfile, a real 1000 Genomes Phase 3 EUR panel (bed/bim/fam, 503 individuals, ~22.7M variants) is enough -- confirmed end-to-end: `ld_clump()` via the real plink1.9 binary against a real 1KG EUR bfile correctly dropped an LD-redundant instrument (r2 = 0.96 with the index SNP, confirmed independently via `plink --r2`) while keeping the index SNP and an independent-locus SNP. A prebuilt bfile is linked from https://mrcieu.github.io/ieugwasr/, or use any local 1000G Phase 3 EUR PLINK build.
