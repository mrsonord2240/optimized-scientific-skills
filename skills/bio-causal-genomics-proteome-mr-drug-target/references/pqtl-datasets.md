# pQTL Datasets and Platform Versioning

Moved verbatim from SKILL.md. Read when choosing or citing a pQTL source, checking whether a target is on a platform, or pinning a release.

## Data Source Taxonomy

| pQTL dataset | Platform | N | Proteins | Reference | Fails when |
|--------------|----------|---|----------|-----------|------------|
| UKB-PPP | Olink Explore 3072 (antibody PEA) | 54,219 | 2923 (54k primary cis-pQTLs across studies; 14,287 primary pQTLs across cis+trans) | Sun 2023 Nature 622:329 | Target not on Olink panel; ancestry mostly EUR; antibody epitope may miss isoforms |
| deCODE | SomaScan v4 (aptamer SOMAmer) | 35,559 | 4907 | Ferkingstad 2021 Nat Genet 53:1712 | Population isolate; LD differs from outbred EUR; aptamers can be PAV-confounded |
| Fenland | SomaScan v4 | 10,708 | 4775 | Pietzner 2021 Science 374:eabj1541 | UK Fenland-specific ascertainment; SomaScan caveats |
| INTERVAL | SomaScan v3 (older) | 3301 | 3622 | Sun 2018 Nature 558:73 | Smaller N; older SomaScan version; useful only as replication |
| ARIC | SomaScan v4 | ~7000 (multi-ancestry sub-cohorts) | 4877 | Zhang 2022 Nat Genet 54:593 | Stratify by ancestry; do not pool |
| FinnGen-PPP | Olink Explore | ~12,000 | ~3000 | FinnGen DF12 (2024 release) | Finnish-specific allele frequencies; not always meta-analyse with UKB |
| AGES-Reykjavik | SomaScan v4 | ~5400 | 4782 | Emilsson 2018 Science 361:769 | Elderly Icelandic cohort; ascertainment bias |
| Olink Explore 1536 disease-specific cohorts (CARDIoGRAMplusC4D, etc) | Olink Explore subset | varies | varies | per-study | Lower N per cohort; use as replication |

Methodology evolves; check the UKB-PPP portal (ukb-ppp.gwas.eu), deCODE Genetics summary-stat releases, and the eQTL Catalogue / GTEx Portal for the current data version. UKB-PPP was substantially re-released in 2024 (extended ancestry meta-analyses); pin a download date in the methods.

### Platform and Cohort Versioning (2024-2026)

- **UKB-PPP Phase 2 (2024)** -- cross-ancestry meta-analyses; ~54k EUR plus multi-ancestry expansion; file layouts differ from Phase 1 (per-protein parquet vs flat TSV); the 14,287 primary pQTL count from Phase 1 shifts in Phase 2 results. Verify the freeze used.
- **FinnGen-PPP** -- Olink Explore platform; DF12 (2024) release is current; Finnish-specific allele frequencies require ancestry-aware downstream analysis.
- **AoU pQTL (2024)** -- All-of-Us proteomics; pre-symptomatic-cohort design strength for reverse-causation control and longitudinal follow-up.
- **SomaScan v4 vs v5** -- v5 expands to ~11k SOMAmers (vs ~5k in v4); binding consistency for shared SOMAmers documented in deCODE / SomaLogic technical notes but not guaranteed; pin platform version.
- **Olink Explore HT** -- ~5,400 proteins (vs ~3,072 in Explore Expansion, ~1,536 in Explore 1536); UKB-PPP Phase 1 used Explore 3072 -- do not assume Explore HT coverage when reading Phase 1 papers.
- Pin specific platform version + release date in the methods section of every cis-MR drug-target manuscript.

