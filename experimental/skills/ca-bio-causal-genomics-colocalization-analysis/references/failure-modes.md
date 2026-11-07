# Per-Method Failure Modes
### coloc.abf -- PP.H3 inflation under multiple causal variants

**Trigger:** Locus has 2+ independent causal signals in moderate LD (r2 ~ 0.3-0.6), one shared with the eQTL and one GWAS-only, **and both datasets are powered to see the signals** (eQTL N >= ~400 at a per-SD effect >= ~0.3, or N ~ 5000 at effect ~ 0.12; GWAS N = 6000). Verified 2026-09-21 by simulation (r2 = 0.52, GWAS causal at SNP A + SNP B, eQTL causal at A, 20 replicates): coloc.abf did not blur to a middling PP.H3; it flipped per replicate between PP.H4 ~ 1 and PP.H3 ~ 1 (about 40% of replicates returned PP.H3 > 0.8 for a truly shared signal), depending on which of the two SNPs led in each dataset. Limited power does **not** trigger it: eQTL N <= 200, or small effects (0.05-0.12 at N <= 1000), gave PP.H1/PP.H0 dominance with PP.H3 <= 0.05. A well-powered single-causal locus at r2 0.5 is also decisive (PP.H4 = 1.0 shared, PP.H3 = 1.0 distinct), so r2 alone is not the trigger.

**Mechanism:** The single-causal-variant assumption forces each dataset onto one SNP. When the lead SNP differs between the GWAS (B) and the eQTL (A) because two nearby signals compete, the per-SNP Bayes factors of the two datasets peak at different SNPs and the posterior goes to H3 (distinct causal variants).

**Symptom:** LocusZoom looks convincing, but `result$summary['PP.H3.abf']` is high while PP.H4 is low, and the answer is unstable under resampling or a different eQTL panel. Because the outcome is near-binary (not intermediate), a decisive PP.H3 at moderate LD in a locus with 2+ GWAS signals is not proof of distinct causal variants: check with coloc.susie. Low PP.H4 with PP.H1 dominating is a power problem, not this failure.

**Fix:** Run coloc.susie (or eCAVIAR or PWCoCo) to allow multiple causal variants. If coloc.susie returns multiple credible sets with one pair showing PP.H4 > 0.75, this is real allelic heterogeneity not failure.

### Failure modes in reference files

- **coloc.susie -- LD reference mismatch** (`estimate_s_rss` lambda > 0.05, spurious credible sets): `references/coloc-susie.md`.
- **Lead-SNP swap and window bias** (PP.H4 changes when the window is re-centred): `references/lead-snp-swap.md`.

### coloc default p12 too liberal for trans-eQTL

**Trigger:** Applying `p12 = 1e-5` (the default) to a trans-eQTL / GWAS pair.

**Mechanism:** The default p12 was calibrated for cis-eQTL where biological proximity makes shared causality reasonable. For trans-eQTL, prior probability of shared causality is much lower; uniform p12 over-favours H4.

**Symptom:** PP.H4 > 0.8 reported, but sensitivity() reveals PP.H4 falls below 0.5 for p12 < 1e-5; replication in independent data fails.

**Fix:** Operational definition: "trans" = >5 Mb from TSS or different chromosome. Default p12=1e-5 over-favors H4 for trans (genome-rare biology). For trans: lower p12 to 5e-6 or 1e-6 AND raise PP.H4 threshold to >= 0.8 (compensate for higher FP risk). Cross-reference Vosa 2021 Nat Genet 53:1300 (eQTLGen trans) for empirical patterns.

### MHC / HLA + chr 8 inversion -- single-causal assumption breaks

**Trigger:** Locus within chr6:25-35 Mb (extended MHC, hg38), or chr8:8.1-11.9 Mb (chr 8 inversion, hg38).

**Mechanism:** The MHC contains classical HLA genes with extreme long-range LD (r2 > 0.5 over many Mb), multiple independent causal haplotypes, and structural variation. The chr 8p23.1 inversion similarly produces long-range LD across megabases of polymorphic inversion alleles. The single-causal-variant assumption is biologically wrong in both regions.

**Symptom:** coloc.abf almost always returns PP.H3 or fragmented PP across H1/H2/H3/H4 even when the underlying biology is well-established (e.g. HLA-DRB1 in autoimmune GWAS).

**Fix (MHC):** Use HLA-imputed classical alleles via SNP2HLA / HIBAG / HLA-TAPAS, then HLA-coloc (Butler-Laporte 2024 medRxiv) -- NOT coloc on SNPs in MHC. OR exclude MHC from genome-wide coloc and report HLA association at the haplotype/allele level. **Fix (chr 8 inversion):** Exclude chr8:8.1-11.9 Mb or pre-condition on inversion genotype before coloc. Never report a single coloc PP.H4 in either region without this caveat.

**Enforce programmatically, do not rely on remembering the coordinates:** call `flag_excluded_region()` (`scripts/flag_excluded_region.R`, already called by `scripts/coloc_abf.R` and `scripts/coloc_susie.R`) on the locus before running coloc.abf/coloc.susie, and abort/redirect rather than compute a naive PP.H4 when it returns non-NA.

### Underpowered eQTL (N < 200)

**Trigger:** Small eQTL discovery (e.g. tissue-specific bulk study, N < 200; per-cell-type sc-eQTL).

**Mechanism:** With low N, varbeta is large; the eQTL's per-SNP Bayes factors are flat; the joint likelihood concentrates on H0 or H1 (GWAS-only).

**Symptom:** PP.H0 or PP.H1 dominates; the eQTL panel shows visible signal but coloc cannot resolve causal vs noise.

**Fix:** Use eQTLGen (N ~ 31k whole-blood) or GTEx v8 (N ~ 70-700 per tissue) where possible. For rare cell types, accept the limitation and report the locus as underpowered rather than claim absence of colocalization.

| eQTL N | Coloc viability | Notes |
|--------|-----------------|-------|
| < 200 | Underpowered | PP.H1 dominant; flag |
| 200-500 | Cis only, modest | Single-tissue cis |
| 500-1000 | Good for cis | Most GTEx v8 tissues |
| >= 1000 | Well-powered | Trans accessible |
| >= 10000 | Meta (eQTLGen) | Cross-tissue / sc |

### Reference QTL panel choice

GTEx v8 (838 donors, 49 tissues, 2020) is the current PredictDB-supported standard. GTEx v10 (released 2024) has limited harmonisation and is not yet PredictDB-default. eQTLGen blood meta-eQTL (N ~ 31k) wins on sample size for blood cis-eQTL discovery, beating any single tissue on power. Always pin version in methods (e.g. "GTEx v8 MASHR-EUR, PredictDB release 2022-01").
