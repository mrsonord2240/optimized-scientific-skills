---
name: bio-clinical-databases-gnomad-frequencies
description: Queries gnomAD v4 (807k samples), v3, v2.1.1, and constraint metrics with grpmax FAF95, bottleneck-group exclusion, LOEUF interpretation, SV/CNV/mtDNA catalogs, and Whiffin max-credible-AF framework. Use when filtering rare variants, applying ACMG BS1/BA1, ranking genes by LoF intolerance, or selecting between v2 (GRCh37) and v4 (GRCh38 + 807k samples).
tool_type: python
primary_tool: requests
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: requests 2.31+, hail 0.2.130+, pandas 2.2+, myvariant 1.0+. Current gnomAD release is **v4.1.1 (March 2026)**; v4.1 fixed the v4.0 AN under-counting issue that inflated rare-variant AF estimates by 5-10%, and v4.1.1 updated constraint metrics and annotation downloads.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- Hail: `hl.version()`; pin to >=0.2.130 for v4 schema

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. The gnomAD browser GraphQL API at `https://gnomad.broadinstitute.org/api` is the supported public endpoint; Hail Tables on Google Cloud Storage at `gs://gcp-public-data--gnomad/` are the supported bulk access.

# gnomAD Frequency Queries and Constraint

**'How rare is this variant in the general population?'** -> Pull allele frequency, grpmax FAF95 (the ACMG-grade frequency), LOEUF gene-level constraint, structural variant catalog, mtDNA frequencies, and the appropriate dataset version per use case.

- Python (single variant): GraphQL via `requests.post('https://gnomad.broadinstitute.org/api', json={'query': ..., 'variables': ...})`
- Python (aggregator): `myvariant.MyVariantInfo().getvariant(hgvs, fields=['gnomad_exome', 'gnomad_genome'])` -- gnomAD 2.1.1 (GRCh37) raw AF only, no FAF95 (myvariant metadata, 2026-09-15)
- Python (bulk): `hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.sites.ht')`

**Data governance:** sending variants derived from patients or research participants to a public API, or processing their exomes in cloud Hail, can disclose them. Do it only under the cohort's consent and institutional approvals, in an approved environment.

## v2.1.1 / v3.1.2 / v4.x: When to Use Which

This is the most consequential decision in any gnomAD query. The releases are **not interchangeable**; choice determines what can and cannot be said about a variant.

| Release | Build | Samples | Use when | Fails when |
|---------|-------|---------|----------|-----------|
| **v2.1.1** | GRCh37 | 125,748 exomes + 15,708 genomes | Constraint metrics needed (LOEUF v2 most-validated); GRCh37 native non-negotiable | GRCh38 native cohort; modern rare-variant FAF95 (use v4) |
| **v3.1.2** | GRCh38 | 76,156 genomes (NO exomes) | Non-coding region rare variants on GRCh38; mtDNA frequencies | Exome variants needed (no exomes); 76k cohort smaller than v4 |
| **v4.0/v4.1.1** | GRCh38 | 730,947 exomes + 76,215 genomes = **807,162 total** | Default for rare-variant filtering, FAF95, and current gene constraint | GRCh37 coordinates (lift over or use v2.1.1); cancer-cohort analysis (no TCGA in v4) |

**Critical caveats:**
- v4 genomes are the SAME 76,215 v3 samples reprocessed against GRCh38 with updated pipelines; not independent.
- ~81% of v2 genomes are also in v3; joint v2+v3 meta-analysis must dedupe at sample ID.
- v4 includes 416,555 UK Biobank exomes under a specific collaboration agreement; check use terms.
- **v4 does NOT include TCGA**, so the `non_cancer` subset is unnecessary; the v4 subset is `non_ukb` (excludes UKB exomes for ancestry rebalancing).
- Liftover v2 (GRCh37 -> GRCh38) is NOT equivalent to v4 native; variant representation differs at ~0.5-1% of sites due to assembly fixes.

## v4 Ancestry Groups: popmax -> grpmax Terminology

v4 ancestry groups: **AFR, AMR, ASJ, EAS, FIN, MID, NFE, SAS, AMI, REMAINING**. The **MID (Middle Eastern) group was new in v4**; previously absorbed into "OTH". The **REMAINING** group (31,256 v4 samples) is individuals who did not cluster with any reference; they contribute to overall AF but not to grpmax.

**Terminology shift:** gnomAD documentation and ACMG-facing narrative uses **grpmax** (genetic ancestry group max) -- replacing the older **popmax** ("population max") term -- to disambiguate genetic ancestry from self-reported race/ethnicity. The public GraphQL schema still exposes legacy field names containing `popmax` (e.g. `faf95.popmax`, `faf95.popmax_population`); these are the grpmax values under the modern terminology. Always check the schema version when writing queries; new browsers may rename these fields.

`grpmax_faf95` is the operational ACMG field. It computes the maximum 95% lower-CI allele frequency, **excluding bottleneck groups** (AMI, ASJ, FIN, REMAINING) because pathogenic founder variants in those groups would otherwise falsely trigger BS1/BA1. MID is included in grpmax but is the smallest non-bottleneck group with highest per-allele variance.

## Filtering Allele Frequency (FAF95): The ACMG-Grade AF

Whiffin 2017 *Genet Med* 19:1151 introduced FAF95 = Poisson lower bound of 95% CI for AF. By construction, AF > FAF95; FAF95 is the conservative frequency for ACMG application.

**Max-credible-AF formula:** `(prevalence x heterogeneity x allelic-contribution) / (penetrance x 2)`. Plug in disease parameters to get the gene-specific BA1 / BS1 threshold; compare against `grpmax_faf95`.

| Code | Threshold | Notes |
|------|-----------|-------|
| **BA1** | AF > 5% in any non-bottleneck group | ClinGen SVI default; VCEPs may override (Hearing Loss VCEP uses 0.5%) |
| **BS1** | AF > gene-specific max-credible-AF | Computed per gene via Whiffin formula |
| **PM2_Supporting** | Absent or ultra-rare in gnomAD | Downgraded from PM2_Moderate in SVI 2020 |

Use `grpmax_faf95`, not raw AF, for BS1/BA1 application; this is the ClinGen-recommended approach.

## Constraint Metrics: pLI, LOEUF, missense Z

Karczewski 2020 *Nature* 581:434 defined LOEUF as the upper bound of the 90% CI of observed/expected pLoF count per gene. LOEUF is **recommended over pLI** because it is continuous and accounts for gene size more rigorously.

| Metric | What | Interpretation |
|--------|------|----------------|
| **LOEUF** | Upper bound of 90% CI of LoF observed/expected ratio | Lower = more LoF-intolerant; **first decile (LOEUF < 0.35 v2; < 0.6 v4) = strongly intolerant** |
| **pLI** | Probability LoF intolerant | Still used; gnomAD team recommends LOEUF for ranking |
| **Missense Z** | Z-score of observed-vs-expected missense | Z > 3.09 = 'constrained' set (~p < 0.001, one-tailed) |
| **Missense O/E** | Observed/expected missense ratio | Continuous form of missense Z |

**Critical version mismatch:**
- v2.1.1 constraint metrics published 2020; v4.1.1 constraint metrics were updated **March 2026** (including all coding loci and chromosomes X/Y).
- The browser API returns constraint for chrX genes on GRCh38 (DMD LOEUF 0.235, checked 2026-09-24). If `gnomad_constraint` comes back null for a gene, fall back to v2.1.1 (`reference_genome: GRCh37`) and say so.
- **LOEUF thresholds are version-specific:** v2 < 0.35 was the prior constrained-gene cutoff. For v4.1.1, gnomAD recommends **LOEUF < 0.45** (about the most constrained 15% of MANE Select transcripts); <0.36 is the v4.1.1 first decile and <0.60 is its first quartile. Compare ranks or percentiles across releases, never raw values alone.

## Subsets: non_cancer, non_neuro, controls

| Release | Subset | Removes | Use when |
|---------|--------|---------|----------|
| v2.1.1 | `non_cancer` | TCGA | Cancer-related variant analysis (avoids circularity) |
| v2.1.1 | `non_neuro` | Psychiatric/neuro cohorts | Neuropsychiatric variant analysis |
| v2.1.1 | `controls` | Cases with known disease (~60k samples) | Disease-association calibration |
| v3.1.2 | `non_v2` | v2 overlapping samples | Independent of v2 |
| v3.1.2 | `controls_and_biobanks` | Disease cases retained, biobanks emphasized | Population-level reference |
| v4 | `non_ukb` | UK Biobank exomes | When EUR-skew of UKB problematic |
| v4 | `non_neuro` | Deprecated | -- |
| v4 | `non_cancer` | Unnecessary (no TCGA in v4) | -- |

SV and CNV catalog details (gnomAD-SV v2/v4 release stats, gnomAD-CNV v4) are in the usage guide.

mtDNA catalog details (Laricchia 2022 *Genome Res* 32:569 frequencies and heteroplasmy thresholds) are in the usage guide.

## VEP Version Pinning

Each gnomAD release pins to a VEP version:
- **v4.1.1 Hail Table downloads use VEP 115**; browser VEP display can remain on an older version for GTEx/pext compatibility
- v2.1.1 uses VEP 85

A variant's consequence prediction can flip between v2 and v4 due to MANE Select adoption and transcript-set updates. Always pin VEP version when reproducing gnomAD annotations.

## Decision Tree by Query Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Single variant AF lookup | GraphQL API (myvariant.info carries gnomAD 2.1.1 only, no FAF95) | Lowest latency; returns full per-ancestry breakdown |
| ACMG BS1/BA1 application | `grpmax_faf95` from v4 | The ClinGen-recommended field |
| Gene-level LoF constraint | LOEUF from the GRCh38 (v4) gene query | Larger sample, more stable |
| Gene-level constraint null in v4 | LOEUF from v2.1.1 (`reference_genome: GRCh37`) | Fallback only when v4 returns no constraint |
| Bulk rare-variant filter (cohort-scale) | Hail Table on GCS | No rate limits; full schema |
| SV frequency | gnomAD-SV v4 (WGS) or gnomAD-CNV v4 (exome) | Choose by data type |
| mtDNA frequency | v3.1 mtDNA release (Laricchia 2022) | Only gnomAD release with mtDNA |
| Cancer-variant analysis | v2.1.1 `non_cancer` subset OR v4 (no TCGA) | Avoid TCGA circularity in v2 |
| Comparison across builds | Use canonical SPDI or CA ID, normalize first | Liftover != native |

## Single Variant Query (GraphQL)

**Goal:** Retrieve exome + genome AF, grpmax, FAF95, and per-ancestry breakdown for one variant.

**Approach:** Hit gnomAD's GraphQL API with explicit dataset version; parse the nested response.

```python
import time
import requests

GNOMAD_API = 'https://gnomad.broadinstitute.org/api'
DATASET_BUILD = {'gnomad_r4': 'GRCh38', 'gnomad_r3': 'GRCh38', 'gnomad_r2_1': 'GRCh37'}

def _post_graphql(query, variables, max_retries=5, base_delay=2.0):
    '''POST to the gnomAD GraphQL API with bounded exponential backoff on HTTP 429.

    The browser API rate-limits sustained querying: HTTP 429 (a plain HTML page, not JSON) was
    observed live after one gene-variant-list query plus a handful of single lookups (checked
    2026-09-15). For more than a few dozen variants, use the Hail Table or sites VCF (see Bulk
    Query via Hail below) instead of looping GraphQL calls.
    '''
    for attempt in range(max_retries):
        r = requests.post(GNOMAD_API, json={'query': query, 'variables': variables}, timeout=30)
        if r.status_code == 429:
            time.sleep(base_delay * (2 ** attempt))
            continue
        r.raise_for_status()
        return r.json()
    r.raise_for_status()  # exhausted retries; surface the last response's error
    return r.json()

def query_variant(chrom, pos, ref, alt, build, dataset='gnomad_r4'):
    '''Query gnomAD GraphQL for variant frequency + grpmax FAF95.

    build: 'GRCh38' or 'GRCh37', the build of the coordinates; checked against the dataset
    (gnomad_r4 / gnomad_r3 = GRCh38, gnomad_r2_1 = GRCh37). GRCh37 ids sent to gnomad_r4 return
    "Variant not found", the same answer as a truly absent variant.
    Returns the variant payload, or None when the variant is not in this dataset; raises on any
    other GraphQL error. Retries with backoff on HTTP 429 (see _post_graphql); when batch-querying,
    pace calls with time.sleep(0.2-0.5) between variants to avoid triggering it.
    '''
    if DATASET_BUILD[dataset] != build:
        raise ValueError(f'{dataset} is {DATASET_BUILD[dataset]} but the coordinates are {build}')
    query = '''
    query VariantById($variantId: String!, $dataset: DatasetId!) {
      variant(variantId: $variantId, dataset: $dataset) {
        variant_id
        rsids
        exome {
          ac
          an
          af
          homozygote_count
          filters
          populations { id ac an }
          faf95 { popmax popmax_population }
        }
        genome {
          ac
          an
          af
          homozygote_count
          filters
          populations { id ac an }
          faf95 { popmax popmax_population }
        }
      }
    }
    '''
    variant_id = f'{chrom}-{pos}-{ref}-{alt}'
    body = _post_graphql(query, {'variantId': variant_id, 'dataset': dataset})
    errors = [e.get('message') for e in body.get('errors') or []]
    if errors and errors != ['Variant not found']:
        raise RuntimeError(f'gnomAD GraphQL error for {variant_id} ({dataset}): {errors}')
    return (body.get('data') or {}).get('variant')


def grpmax_faf95(payload):
    '''Extract the grpmax FAF95; the ACMG-grade frequency. Excludes bottleneck groups.

    faf95 is None both for absent variants (source='absent') and for present variants whose FAF95
    is undefined because too few alleles were seen (source='present_faf95_undefined', e.g. AC=1).
    '''
    if payload is None:
        return {'faf95': None, 'grpmax_ancestry': None, 'source': 'absent'}
    for source in ('exome', 'genome'):
        faf = (payload.get(source) or {}).get('faf95') or {}
        if faf.get('popmax') is not None:
            return {'faf95': faf['popmax'], 'grpmax_ancestry': faf.get('popmax_population'), 'source': source}
    return {'faf95': None, 'grpmax_ancestry': None, 'source': 'present_faf95_undefined'}
```

## ACMG BS1/BA1 Application

**Goal:** Apply Whiffin max-credible-AF framework to a candidate variant.

**Approach:** Compute the gene-specific BS1 threshold from disease parameters, compare to `grpmax_faf95`.

```python
def max_credible_af(prevalence, max_allelic_contribution=1.0, max_genetic_contribution=1.0,
                    penetrance=1.0):
    '''Whiffin 2017 max-credible-AF formula.

    Args:
        prevalence: disease prevalence (e.g., 1/10000 = 1e-4)
        max_allelic_contribution: max contribution of single allele to disease in any case
        max_genetic_contribution: max contribution of this gene to disease in any case
        penetrance: probability that variant carriers develop disease

    Returns: max-credible per-allele frequency under dominant inheritance (use /2 for AR)
    '''
    return (prevalence * max_genetic_contribution * max_allelic_contribution) / (penetrance * 2)


def apply_bs1_ba1(grpmax_faf95_val, max_credible, ba1_threshold=0.05):
    '''Apply ClinGen SVI BS1/BA1 criteria.

    BA1 default 5% per ClinGen SVI; VCEP-specific overrides exist (Hearing Loss = 0.5%).
    BS1 = max-credible-AF specific to gene+disease.
    grpmax_faf95_val None = absent, or present with FAF95 undefined (see grpmax_faf95()['source']).
    These are research annotation tags for a variant, not a classification.
    '''
    if grpmax_faf95_val is None:
        return 'PM2_Supporting'  # Absent or ultra-rare
    if grpmax_faf95_val > ba1_threshold:
        return 'BA1'
    if grpmax_faf95_val > max_credible:
        return 'BS1'
    return None  # No criterion triggered; variant is consistent with rare-disease causation
```

## Gene-Level Constraint (LOEUF)

**Goal:** Retrieve gene constraint metrics with awareness of v2/v4 version differences.

**Approach:** Use the GRCh38 (v4) LOEUF; fall back to v2.1.1 only when v4 returns no constraint. Report LOEUF decile, not raw value, to avoid cross-version comparison errors.

```python
def query_gene_constraint(gene_symbol, dataset='gnomad_r4'):
    '''Pull gene constraint metrics (chrX genes included: DMD LOEUF 0.235 on 2026-09-15).'''
    query = '''
    query GeneById($symbol: String!) {
      gene(gene_symbol: $symbol, reference_genome: GRCh38) {
        gene_id
        symbol
        chrom
        gnomad_constraint {
          oe_lof
          oe_lof_lower
          oe_lof_upper
          oe_mis
          oe_mis_upper
          pli
          mis_z
        }
      }
    }
    '''
    gene = _post_graphql(query, {'symbol': gene_symbol}).get('data', {}).get('gene')
    if gene is None:
        return None
    if gene.get('gnomad_constraint') is None:
        gene['constraint_note'] = ('no constraint returned for this gene; try the v2.1.1 values '
                                   'with reference_genome: GRCh37')
    return gene
```

## Bulk Query via Hail (cohort-scale)

**Goal:** Filter millions of variants by AF, grpmax, or LOEUF without API rate limits.

**Approach:** Read gnomAD v4 Hail Table from Google Cloud Storage; use `hl.read_table()` + filter operations.

```python
import hail as hl

def init_hail_for_gnomad():
    '''Initialize Hail for gnomAD v4 GCS access. Requires Hail 0.2.130+.'''
    hl.init(default_reference='GRCh38')


def filter_rare_variants_hail(input_vcf, max_grpmax_faf95=0.0001, output_path='filtered.mt'):
    '''Filter input MT to variants below grpmax FAF95 threshold using gnomAD v4 exomes.'''
    ht_v4 = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/'
                          'gnomad.exomes.v4.1.sites.ht')
    mt = hl.import_vcf(input_vcf, reference_genome='GRCh38')
    mt = mt.annotate_rows(gnomad=ht_v4[mt.locus, mt.alleles])
    mt = mt.filter_rows(
        (hl.is_missing(mt.gnomad.grpmax_faf95)) |
        (mt.gnomad.grpmax_faf95.faf95 < max_grpmax_faf95)
    )
    mt.write(output_path, overwrite=True)
    return mt
```

## Per-Operation Failure Modes

**1. Using popmax/AF where grpmax_faf95 belongs**
- Trigger: Apply BS1 with raw AF instead of FAF95.
- Mechanism: Raw AF inflates for low-N populations; FAF95 is the lower-bound CI; conservative.
- Symptom: Pathogenic variants falsely categorized BS1 in small-N ancestry groups (especially MID with v4's smallest sample size).
- Fix: Use `grpmax_faf95.popmax` field; not `populations[i].af`.

**2. Failing to exclude bottleneck groups**
- Trigger: Compute grpmax including AMI, ASJ, FIN, REMAINING.
- Mechanism: Founder variants in bottleneck groups can reach AF > 5% but are not population-general; would falsely trigger BA1.
- Symptom: Founder-population pathogenic variants reported benign.
- Fix: Use gnomAD's pre-computed `grpmax_faf95` which excludes bottleneck groups by design.

**3. GRCh37 coordinates sent to a GRCh38 dataset**
- Trigger: Query `gnomad_r4` with GRCh37 positions (e.g. rs334 as 11-5248232-T-A).
- Mechanism: The API answers `{"errors": [{"message": "Variant not found"}], "data": {"variant": null}}`; reading only `data` turns a build mismatch into "absent".
- Symptom: A common variant is reported absent from gnomAD.
- Fix: Check the build first; surface the GraphQL `errors` array; query `gnomad_r2_1` for GRCh37 or lift over to GRCh38 (rs334 = 11-5227002-T-A).
- Residual case: `query_variant`'s `build` check only catches a caller who contradicts themselves (`build='GRCh37'` with `dataset='gnomad_r4'`); it cannot catch coordinates that are simply wrong for a `build` the caller states correctly. rs334's GRCh37 coordinates (11-5248232-T-A) queried against `gnomad_r4` with `build='GRCh38'` still return the identical `Variant not found` (confirmed live 2026-09-16) -- the API cannot tell a wrong-build id from a truly absent one. Before trusting a batch of "absent" results, spot-check one known common variant (rs334 -> GRCh38 11-5227002-T-A, genome AF ~1.3%) or confirm the reference allele via NCBI Variation Services `/v0/spdi/{seq_id}:{pos-1}:{ref}:{alt}/canonical_representative` (SPDI positions are 0-based; a `Disambiguation exception` warning there means the asserted reference doesn't match the named assembly -- a build-mismatch signal) or `/v0/refsnp/{rsid}`.

**4. Comparing LOEUF absolute values across v2/v4**
- Trigger: "v4 LOEUF for GENE-X is 0.45; v2 was 0.30; has it become more tolerant?"
- Mechanism: Larger v4 sample shifts the LOEUF distribution upward; first-decile threshold shifted v2 < 0.35 -> v4 < 0.6.
- Symptom: Genes appear to lose constraint between versions when they have not.
- Fix: Compare deciles, not absolute values; or stay within one version.

**5. v2 -> v4 liftover assumed equivalent**
- Trigger: Project v2 GRCh37 variants onto GRCh38 with CrossMap, treat as v4 native.
- Mechanism: ~0.5-1% of sites have different representations after liftover due to assembly fixes (e.g., gaps closed, contigs joined).
- Symptom: Inconsistent AFs at low rate; failed cross-version reproducibility.
- Fix: Query v4 native by GRCh38 coordinates directly; do not use liftover output as v4-equivalent.

**6. UKB sample contamination of grpmax**
- Trigger: Compute grpmax across v4 default subset; observe inflated NFE/SAS.
- Mechanism: 416,555 UK Biobank exomes dominate the v4 NFE+SAS subsets.
- Symptom: Variants common in UKB but rare globally falsely look common.
- Fix: Use `non_ukb` subset for grpmax when ancestry composition matters.

**7. v3 vs v4 confusion; "I want WGS"**
- Trigger: User says "I want WGS AFs" and pipeline pulls v4 genomes.
- Mechanism: v4 genomes are the SAME 76,215 v3 samples reprocessed against GRCh38; not independent.
- Symptom: WGS AFs appear identical to v3.1.2; not a bug, but worth flagging.
- Fix: Document that v4 genomes = v3 genomes reprocessed; for true independent WGS, no such resource yet exists at scale.

**8. Constraint applied to multi-isoform gene without transcript awareness**
- Trigger: Apply LOEUF "for the gene" when LoF is isoform-specific.
- Mechanism: gnomAD constraint is computed on the canonical transcript; tissue-specific or alternative isoforms may have different LoF tolerance.
- Symptom: Mis-prioritization of variants on minor transcripts.
- Fix: Cross-check with MANE Select; for isoform-specific LoF, use isoform-level constraint where available (rare).

## Reconciliation: When Sources Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| ClinVar P vs gnomAD `grpmax_faf95` > 1% | Founder-population pathogenic; or ClinVar is stale low-star | Apply Whiffin max-credible-AF for the gene; check ClinVar star/freshness |
| v2 LOEUF < 0.35 vs v4 LOEUF = 0.5 | Distribution shifted with v4 sample size, not biology | Use deciles; v4 first decile = < 0.6 |
| v2 AF != v4 AF for same variant | Sample overlap (v3 in v4) + new exomes; expected | Trust v4 default; non-overlapping subsets via `non_v2` or `non_ukb` |
| Variant present in v3 genomes, absent v4 exomes | Variant outside exome capture region (intronic, intergenic) | Use v3.1.2 or v4 genomes for non-coding |
| gnomAD-SV v2 vs v4 different breakpoints | v2 GRCh37, v4 GRCh38; assembly fixes shift coords | Use v4 native; document build |
| Browser shows lower AF than Hail Table | Browser pre-filters with `filters=PASS`; Hail Table includes all | Apply `filters` filter in Hail explicitly |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| BA1 default | grpmax_faf95 > 5% in non-bottleneck group | Richards 2015 + ClinGen SVI |
| BS1 | grpmax_faf95 > gene-specific max-credible-AF | Whiffin 2017 |
| PM2_Supporting | Absent or ultra-rare in gnomAD | SVI 2020 downgrade |
| LOEUF first decile v2 | < 0.35 | Karczewski 2020 |
| LOEUF first decile v4.1.1 | < 0.36 | gnomAD v4.1.1 release March 2026 |
| Recommended constrained-gene cutoff v4.1.1 | < 0.45 | gnomAD v4.1.1 release March 2026 |
| Missense Z constrained | Z > 3.09 (~p < 0.001) | Samocha 2014 |
| mtDNA heteroplasmy carrier threshold | >=10% heteroplasmy | Laricchia 2022 |
| v4 sample size | 730,947 exomes + 76,215 genomes = 807,162 | gnomAD v4.0 release Nov 2023 |
| Bottleneck groups (excluded from grpmax) | AMI, ASJ, FIN, REMAINING | gnomAD v4 documentation |
| API rate limit | None published; ~10 req/s practical | gnomAD browser GraphQL |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| `Cannot read property 'af' of undefined` | Variant not in dataset; `variant` returned null | Check `if payload is None`; absence is biologically informative |
| FAF95 is null for a present variant | FAF95 is undefined when too few alleles are observed (AC=1 exome: `popmax: null`) | Check AC and AN directly; do not read null as absent |
| `variant` null with `errors: Variant not found` | Not in this dataset, or GRCh37 coordinates sent to a GRCh38 dataset | Confirm the build; `query_variant` requires it |
| Variant filter status `AC0` or `RF` | Failed gnomAD QC | Variants with non-`PASS` should usually be excluded from analysis |
| Different AFs between gnomAD browser and Hail Table | Browser auto-applies PASS filter; Hail does not | Filter `filters.size() == 0` (i.e., `PASS`) in Hail |
| LOEUF appears worse in v4 vs v2 | Distribution shifted with larger sample | Compare deciles, not absolute values |
| SV not found in v4-SV | v2-SV is GRCh37, v4-SV is GRCh38; or variant not called in WGS | Try v2-SV with liftover; or check gnomAD-CNV for exome-derived |
| mtDNA variant missing | Only v3.1 has mtDNA; not in v4 | Query v3.1 directly |

Anticipated reviewer pushback and standard responses are in the usage guide.

## References

- Chen S et al. 2024. A genomic mutational constraint map using variation in 76,156 human genomes. *Nature* 625:92.
- Karczewski KJ et al. 2020. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 581:434.
- Samocha KE et al. 2014. A framework for the interpretation of de novo mutation in human disease. *Nat Genet* 46:944.
- Collins RL et al. 2020. A structural variation reference for medical and population genetics. *Nature* 581:444.
- Laricchia KM et al. 2022. Mitochondrial DNA variation across 56,434 individuals in gnomAD. *Genome Res* 32:569.
- Whiffin N et al. 2017. Using high-resolution variant frequencies to empower clinical genome interpretation. *Genet Med* 19:1151.
- ClinGen guidance on gnomAD v4 (March 2024): `https://clinicalgenome.org/site/assets/files/9445/clingen_guidance_to_vceps_regarding_the_use_of_gnomad_v4_march_2024.pdf`
- gnomAD v4 release notes: `https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/`
- gnomAD v4.1.1 constraint update: `https://gnomad.broadinstitute.org/news/2026-03-gnomad-v4-1-1/`

## Related Skills

- clinical-databases/clinvar-lookup - Pathogenicity classification (gnomAD AF used for BS1/BA1)
- clinical-databases/acmg-classification - Whiffin FAF95 framework applied to ACMG criteria
- clinical-databases/variant-prioritization - Rare-disease pipeline using grpmax_faf95
- clinical-databases/myvariant-queries - Aggregated queries including gnomAD overlay
- population-genetics/population-structure - Population stratification background
