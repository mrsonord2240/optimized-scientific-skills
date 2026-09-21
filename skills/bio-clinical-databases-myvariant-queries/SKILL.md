---
name: bio-clinical-databases-myvariant-queries
description: Queries myvariant.info BioThings aggregator for ClinVar, gnomAD, dbSNP, dbNSFP, COSMIC, CADD, and CIViC annotations in batched, version-tracked requests. Use when annotating variant lists from multiple databases simultaneously without managing per-source APIs, and when reproducibility-grade analyses require recording source data versions from /v1/metadata.
tool_type: python
primary_tool: myvariant
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: myvariant 1.0.0+, requests 2.31+, pandas 2.2+. myvariant.info aggregates >=21 sources; the operative version of each source is queryable via the `/v1/metadata` endpoint (records carry no `_meta`). Field paths and versions below were checked against the live service on 2026-09-15.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. dbNSFP version drift is the dominant staleness vector: AlphaMissense was added to dbNSFP v4.4 (~2024); querying `dbnsfp.alphamissense.score` returns whatever version of dbNSFP is currently loaded; check `src.dbnsfp.version` in `/v1/metadata` (4.8a on 2026-09-15).

# MyVariant.info Queries; Aggregated Annotation

**'Annotate my variants with ClinVar + gnomAD + CADD + AlphaMissense in one batch'** -> Query the BioThings myvariant.info aggregator with field selection and version tracking, then parse nested responses.

- Python: `myvariant.MyVariantInfo().getvariant(hgvs_or_rsid, fields=['clinvar', 'gnomad_exome', 'dbnsfp'])`
- Python (batch): `mv.getvariants(ids_list, fields=...)`; the client chunks the list itself (a raw POST over 1,000 IDs returns HTTP 400)
- Python (search): `mv.query('clinvar.gene.symbol:BRCA1 AND clinvar.rcv.clinical_significance:Pathogenic')`
- REST: `GET https://myvariant.info/v1/variant/{hgvs_or_id}?fields=...`
- Bulk: `POST https://myvariant.info/v1/variant` with comma-separated IDs

## BioThings Architecture (Lelong 2022 *Bioinformatics*)

myvariant.info is one of three flagship BioThings APIs (with MyGene.info and MyChem.info). All three share the BioThings SDK, which auto-deploys an Elasticsearch index from heterogeneous source files via per-source dataloaders. The 2022 paper formalized the SDK; the architecture itself is older (Xin 2016 *Genome Biol*).

- Elasticsearch-backed: queries use Lucene operators (AND, OR, NOT, range like `cadd.phred:>20`)
- Dotted-field-name syntax for nested JSON
- The `_id` field is canonical HGVS-g per record (e.g., `chr7:g.117199644G>A`)
- `_id` coordinates are **hg19/GRCh37**, not GRCh38. A GRCh38 HGVS-g string returns 404/notfound by `_id` even when the variant is loaded (live-confirmed 2026-09-16: `chr17:g.43106487A>C` (GRCh38) -> 404; the same BRCA1 variant is keyed `chr17:g.41258504A>C` (hg19)). For GRCh38 input, either convert to hg19 first (rsID via `clinical-databases/dbsnp-queries`, or SPDI via NCBI Variation Services / `clinical-databases/clinvar-lookup`), or search the hg38-namespaced fields directly: `clinvar.hg38.start`/`clinvar.hg38.end` or `dbnsfp.hg38.start`/`dbnsfp.hg38.end` (live-confirmed working 2026-09-16 — do not use the unprefixed `hg38.start`, which is not a valid field; see `find_high_cadd_in_region` below, where CADD has no hg38 coordinates at all and hg19 is the only option)

## Aggregated Sources: ~21 and Counting

| Source | What | Notes |
|--------|------|-------|
| ClinVar | Pathogenicity under `clinvar.rcv[]` | Loaded release 2025-05 (not weekly); current ClinVar -> clinvar-lookup |
| gnomAD 2.1.1 exomes + genomes (GRCh37) | Population AF (`gnomad_exome.af.af`) | No FAF95 fields; v4 grpmax FAF95 -> gnomad-frequencies |
| dbSNP Build 156 | rsID + alleles | RsMergeArch resolved |
| dbNSFP 4.8a | Meta-aggregator of 40+ in silico predictors | Includes AlphaMissense, REVEL, BayesDel |
| CADD | Deleteriousness | Genome-wide; field `cadd.phred` (not `dbnsfp.cadd`) |
| CIViC | Cancer interpretation | Per-disease |
| COSMIC | Somatic variants | Catalogue of Somatic Mutations |
| EVS | Exome Variant Server | Legacy (deprecated by gnomAD) |
| ExAC | ExAC frequencies | Legacy (superseded by gnomAD) |
| GRASP | GWAS associations | -- |
| GWAS Catalog | Curated GWAS | -- |
| Wellderly | Disease-resistant elderly cohort | -- |
| EMV | -- | -- |
| DOCM | Database of Curated Mutations | -- |
| ICGC | International cancer | -- |
| MutDB | -- | -- |
| GO | Gene Ontology | -- |
| Snpeff | snpEff annotations | -- |
| GeneReviews | Disease/gene reviews | -- |
| MutPred | Functional impact | -- |

**dbNSFP is itself an aggregator.** Querying `dbnsfp.alphamissense.score` returns the version that dbNSFP loaded, not AlphaMissense direct. The lag from publication (Cheng 2023 *Science*) to integration into myvariant.info is typically 6-18 months via dbNSFP.

## Scopes and Query Forms

| Endpoint | Method | Use |
|----------|--------|-----|
| `/v1/variant/{id}` | GET | Single canonical-ID lookup |
| `/v1/variant` | POST (batched IDs) | Batch lookup, max 1,000 IDs per POST (HTTP 400 above) |
| `/v1/query?q={lucene}` | GET | Flexible Elasticsearch search |
| `/v1/metadata` | GET | Per-source versions |

**Scopes** (the `scopes` parameter on `/v1/query` POST) specifies which fields to match an input ID against: `hgvs`, `rsid`, `dbsnp.rsid`, `dbnsfp.genename`, `chrom`, `_id`. The `_id` is canonical HGVS-g.

## Reproducibility: Source Versions from `/v1/metadata`

Records carry no per-record `_meta`; read the instance-wide source versions once per run:

```python
import requests

meta = requests.get('https://myvariant.info/v1/metadata', timeout=30).json()
versions = {name: src.get('version') for name, src in meta.get('src', {}).items()}
print(versions['dbnsfp'], versions['clinvar'], versions['gnomad'])  # '4.8a' '2025-05' '2.1.1' on 2026-09-15
```

For reproducibility, record per-source versions in analysis output alongside results.

## Comparison to Alternatives

| Tool | Approach | When to use |
|------|----------|-------------|
| **myvariant.info** | Cloud aggregator, ES-backed | Quick batch annotation, no local setup |
| **OpenCRAVAT** (Pagel 2020 *JCO Clin Cancer Inform*) | Local install, modular annotators | Offline / PHI-sensitive |
| **VarSome** (Kopanos 2019 *Bioinformatics* 35:1978; commercial) | Hosted, 22 sources | 82% ACMG criteria auto-application (highest); clinical labs |
| **Franklin / Genoox** | Commercial hosted | 59 data sources; family/cohort analysis |
| **GeneBe.net** (Stawiński 2024 *Clin Genet*) | Open-source web + API | Free Tavtigian-point-system-based ACMG; comparable to VarSome |
| **ANNOVAR / VEP / snpEff** | Local annotation tools | Pipeline integration, batch annotation, no ACMG |

**myvariant.info does NOT produce ACMG calls**; it is purely an annotation aggregator. Pair with InterVar, GeneBe, or the `acmg-classification` skill for classification.

## Decision Tree by Query Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Single variant batch annotation | `getvariant(hgvs, fields=...)` | One call, all aggregated sources |
| 10-1000 variants | `getvariants(list, fields=...)` | Batch endpoint, up to 1000 |
| > 1000 variants | `getvariants(list, fields=...)` (client chunks); raw POST must chunk to 1000 | A POST over 1,000 IDs returns HTTP 400 |
| Search by gene + pathogenicity | `mv.query('clinvar.gene.symbol:BRCA1 AND clinvar.rcv.clinical_significance:Pathogenic', size=200)`; read `total` | Elasticsearch Lucene; `hits` holds at most `size` records |
| ACMG-grade pipeline | myvariant for annotation -> InterVar / GeneBe for classification | myvariant does not produce ACMG calls |
| Offline / PHI-sensitive | OpenCRAVAT or VEP locally | myvariant requires HTTP |
| Reproducibility | Record `/v1/metadata` source versions | dbNSFP version is the dominant staleness vector |
| Source-specific deep dive | Use the source-specific skill (clinvar-lookup, gnomad-frequencies) | myvariant is aggregator-grade, not source-deep |

## Standard Annotation Workflow

**Goal:** Annotate a list of variants with the canonical clinical fields for downstream prioritization.

**Approach:** Batch `getvariants` with explicit field list; record `/v1/metadata` versions; convert to DataFrame.

**Data governance:** batch-annotating variant lists derived from patients or research participants sends them to myvariant.info, a public API. Confirm consent and IRB/data-use approvals cover this before batch-querying; route PHI-sensitive work to OpenCRAVAT or another local tool instead (see Comparison to Alternatives) when those approvals aren't in place.

```python
import myvariant
import pandas as pd
import requests

mv = myvariant.MyVariantInfo()

CLINICAL_FIELDS = [
    'clinvar.rcv.clinical_significance',   # ClinVar significance lives under clinvar.rcv[]
    'clinvar.rcv.review_status',
    'clinvar.variant_id',
    'gnomad_exome.af.af',                  # gnomAD 2.1.1 (GRCh37); myvariant has no FAF95
    'gnomad_exome.an.an',
    'gnomad_genome.af.af',
    'dbsnp.rsid',
    'dbnsfp.alphamissense.score',
    'dbnsfp.alphamissense.pred',
    'dbnsfp.revel.score',
    'cadd.phred',                          # top-level CADD source, not dbnsfp.cadd
    'cosmic.cosmic_id',                    # no SpliceAI field in myvariant (only cadd.dst2splice)
    'civic.openCravatUrl'
]

def source_versions():
    '''Instance-wide per-source versions; records carry no _meta.'''
    meta = requests.get('https://myvariant.info/v1/metadata', timeout=30).json()
    return {name: src.get('version') for name, src in meta.get('src', {}).items()}

def _as_list(value):
    return value if isinstance(value, list) else ([value] if value is not None else [])

def annotate_variant_list(hgvs_list):
    '''Batch-annotate variants with ClinVar / gnomAD 2.1.1 / dbNSFP / CADD / COSMIC fields.

    getvariants() chunks the list itself. A multi-allelic rsID returns one record per allele.
    '''
    rows = []
    for r in mv.getvariants(hgvs_list, fields=CLINICAL_FIELDS):
        if r.get('notfound'):
            rows.append({'variant': r.get('query'), 'notfound': True})
            continue
        rcv = _as_list((r.get('clinvar') or {}).get('rcv'))
        gnomad_e = r.get('gnomad_exome') or {}
        gnomad_g = r.get('gnomad_genome') or {}
        dbnsfp = r.get('dbnsfp') or {}
        rows.append({
            'variant': r.get('query'),
            'myvariant_id': r.get('_id'),
            'clinvar_sig': sorted({x.get('clinical_significance') for x in rcv if x.get('clinical_significance')}),
            'clinvar_review': sorted({x.get('review_status') for x in rcv if x.get('review_status')}),
            'gnomad_v2_af': (gnomad_e.get('af') or {}).get('af') or (gnomad_g.get('af') or {}).get('af'),
            'rsid': (r.get('dbsnp') or {}).get('rsid'),
            'alphamissense': (dbnsfp.get('alphamissense') or {}).get('score'),
            'revel': (dbnsfp.get('revel') or {}).get('score'),
            'cadd_phred': (r.get('cadd') or {}).get('phred'),
            'cosmic_id': (r.get('cosmic') or {}).get('cosmic_id')
        })
    return pd.DataFrame(rows), source_versions()
```

## Elasticsearch Query Patterns

**Goal:** Search beyond canonical IDs; e.g., all pathogenic variants in a gene, all variants in a genomic region with CADD > 20.

**Approach:** Lucene syntax in `mv.query()`; support boolean operators, ranges, wildcards.

```python
def find_pathogenic_in_gene(gene_symbol, max_results=500):
    '''Find ClinVar P/LP variants in a gene. Returns (total, hits); hits holds at most max_results
    records (BRCA1: total 3,933 on 2026-09-15), so compare len(hits) with total.'''
    query = f'clinvar.gene.symbol:{gene_symbol} AND '\
            'clinvar.rcv.clinical_significance:(Pathogenic OR "Likely pathogenic")'
    res = mv.query(query, size=max_results, fields=['_id', 'clinvar.rcv.clinical_significance',
                                                      'clinvar.rcv.review_status'])
    return res.get('total'), res.get('hits', [])


def find_high_cadd_in_region(chrom, start, end, min_cadd=25):
    '''Find variants in region with CADD phred above threshold.

    start/end are GRCh37 (hg19) positions. hg38.start combined with cadd.phred returned 0 hits
    for the BRAF V600E window on 2026-09-15 (145 with hg19.start), so query CADD by hg19 coordinates.
    '''
    query = f'chrom:{chrom} AND hg19.start:[{start} TO {end}] AND '\
            f'cadd.phred:>{min_cadd}'
    return mv.query(query, size=500, fields=['_id', 'cadd.phred', 'clinvar.rcv.clinical_significance'])


def find_alphamissense_pathogenic(gene, min_score=0.564):
    '''Find AlphaMissense pathogenic missense in a gene.

    Note: Cheng 2023 dev cutoff is 0.564 BUT this is NOT the Pejaver-style calibrated
    PP3 threshold. ClinGen has not endorsed AlphaMissense thresholds as of May 2026;
    use AlphaMissense as supporting evidence only.
    '''
    query = f'dbnsfp.genename:{gene} AND dbnsfp.alphamissense.score:>{min_score}'
    return mv.query(query, size=500, fields=['_id', 'dbnsfp.alphamissense', 'clinvar.rcv.clinical_significance'])
```

## Per-Operation Failure Modes

**1. Stale dbNSFP version**
- Trigger: Query `dbnsfp.alphamissense.score` and trust as current.
- Mechanism: dbNSFP version lags new tools by 6-18 months; AlphaMissense (Cheng 2023) was integrated into dbNSFP 4.4 (~2024).
- Symptom: Predictions appear missing for variants newly scored by AlphaMissense.
- Fix: Check `src.dbnsfp.version` in `/v1/metadata`; for cutting-edge predictions query AlphaMissense API directly.

**2. Treating AlphaMissense dev threshold as PP3-calibrated**
- Trigger: Apply AlphaMissense score >0.564 as PP3 evidence.
- Mechanism: Cheng 2023 developer-recommended threshold is NOT the Pejaver 2022 calibration framework; ClinGen has not endorsed strength-graded thresholds for AlphaMissense as of May 2026.
- Symptom: Over-application of PP3 in ACMG classification.
- Fix: Treat AlphaMissense as supporting evidence; defer to `clinical-databases/acmg-classification` for calibrated thresholds.

**3. Stacking REVEL + BayesDel + AlphaMissense as independent evidence**
- Trigger: Use multiple in silico predictors as additive PP3 evidence.
- Mechanism: REVEL, BayesDel, VEST4 share ClinVar/HGMD training labels; AlphaMissense is partially independent but correlates strongly with conservation.
- Symptom: Inflated PP3 strength; double-counting.
- Fix: Apply ONE predictor per variant (Pejaver 2022 explicit recommendation).

**4. Rate-limit ignorance**
- Trigger: Loop sequentially over 10k variants with `getvariant()`.
- Mechanism: ~1000 req/sec aggregate per source IP; loops trigger throttling or 429s.
- Symptom: Increasing latency, eventual failures.
- Fix: Use `getvariants(chunk, fields=...)` with chunk size 1000; sleep ~0.5s between chunks.

**5. Field-path errors silently return None**
- Trigger: Query a path that is not a myvariant field, e.g. `clinvar.clinical_significance`, `dbnsfp.cadd.phred` or `gnomad_exome.faf95`.
- Mechanism: Elasticsearch returns null for non-existent paths and 0 hits for searches on them; no error raised.
- Symptom: All values None, or a search total of 0 (BRCA1 P/LP: 0 with `clinvar.clinical_significance`, 3,933 with `clinvar.rcv.clinical_significance`).
- Fix: Check `/v1/metadata/fields`; ClinVar significance and review status are under `clinvar.rcv`, CADD is `cadd.phred`; confirm a known variant (BRAF V600E) returns non-empty values.

**6. Multi-allelic rsID returns one variant only**
- Trigger: Query `rs12345` and treat returned variant as the variant of interest.
- Mechanism: rsID is a cluster identifier; multi-allelic clusters return multiple records.
- Symptom: Wrong allele returned for ~6-8% of rsIDs.
- Fix: Use HGVS-g instead of rsID for unambiguous lookups; for rsID queries, inspect all returned variants and filter by allele.

**7. Sample overlap between sources**
- Trigger: Treat ClinVar + gnomAD as independent corroboration.
- Mechanism: ClinVar variants often *derive from* gnomAD-included individuals; not statistically independent.
- Symptom: False sense of independent evidence in ACMG application.
- Fix: ClinVar P + gnomAD rare is operationally additive evidence per ACMG but understand the populations may overlap.

## Reconciliation: When Sources Inside myvariant Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| dbNSFP REVEL != ClinVar PP3 strength | Different curation cohort | Use Pejaver 2022 calibrated thresholds (see `acmg-classification`) |
| ClinVar P + AlphaMissense benign | NMD-escape region, alternative isoform, ClinVar P stale | Cross-check with conservation, splicing predictions |
| gnomAD AF differs across exome vs genome | Sample sizes differ; myvariant carries gnomAD 2.1.1 (125,748 exomes, 15,708 genomes) | Prefer exome AF; use gnomad-frequencies for v4 FAF95 |
| COSMIC + ClinVar overlap | True dual-classification (germline + somatic) | Report both contexts |
| Variant missing from one source | Source-specific coverage gaps | Cross-check directly with primary source skill (clinvar-lookup, gnomad-frequencies) |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| Batch endpoint cap | 1000 IDs per POST; HTTP 400 above (the Python client chunks) | live check 2026-09-15 |
| Rate limit | ~1000 req/sec aggregate; lower per IP | myvariant.info docs |
| dbNSFP refresh lag | 6-18 months from primary source release | dbNSFP release history |
| `/v1/metadata` | Per-source versions for the instance; records carry no `_meta` | live check 2026-09-15 |
| Chrom:pos term | Quote a full HGVS id for an exact match (`"chr7:g.140453136A>T"`, 1 hit); the bare `chrom:pos` term also works unescaped (`chr7:140453136`, 5 hits, all variants at that position) | live check 2026-09-16 |
| Multi-allelic rsID | ~6-8% of dbSNP rsIDs are multi-allelic | operational estimate |
| AlphaMissense PP3 calibration | NOT yet ClinGen-endorsed (as of May 2026) | ClinGen SVI |
| REVEL PP3_Strong calibration | >= 0.932 per Pejaver 2022 | Pejaver 2022 *AJHG* |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| `KeyError: 'gnomad_exome'` | Variant absent from gnomAD exome dataset | Use `.get('gnomad_exome', {})` defensively |
| `None` for AlphaMissense on rare variants | dbNSFP coverage gap; variant in alt-spliced isoform | Query AlphaMissense API directly, or accept None |
| Search returns 0 hits despite known matches | `:` in a `chrom:pos` term was escaped (`chr7\:140453136`); myvariant does not want that escape | Use the unescaped `chrom:pos` term, or quote a full HGVS id for one exact variant |
| Batch returns < input IDs | Some IDs not in any source | Check `notfound` field in response |
| Different AF in myvariant vs gnomAD browser | myvariant carries gnomAD 2.1.1; the browser defaults to v4 | Check `src.gnomad.version` in `/v1/metadata`; use gnomad-frequencies for v4 |
| Search returns 0 hits for ClinVar or CADD terms | Non-existent field path | Use `clinvar.rcv.clinical_significance`, `clinvar.rcv.review_status`, `cadd.phred` |
| 503 on bulk query | Rate limit | Reduce chunk to 500; sleep 1s between |
| `_id` doesn't match input | myvariant uses canonical HGVS-g; input was rsID or non-canonical | Re-query by `_id` after first resolution |
| GRCh38 HGVS-g by `_id` returns 404/notfound | `_id` is hg19/GRCh37; myvariant does not index GRCh38 coordinates as `_id` | Query `clinvar.hg38.start`/`dbnsfp.hg38.start` to find the hg19 `_id`, or convert via rsID/SPDI first |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "myvariant.info is just an aggregator; why not query sources directly?" | Aggregator avoids per-source API setup; sufficient for single-source unique queries we defer to source-specific skills. |
| "This annotation differs from VarSome" | VarSome uses its own ACMG implementation; myvariant.info does NOT produce ACMG calls; we pair with `acmg-classification`. |
| "dbNSFP REVEL differs from REVEL website" | dbNSFP version is on the order of 1 year behind primary; check `src.dbnsfp.version` in `/v1/metadata`. |
| "AlphaMissense calibration thresholds were missed" | AlphaMissense is integrated via dbNSFP; PP3 calibration is in `clinical-databases/acmg-classification` skill. |
| "Why not OpenCRAVAT?" | OpenCRAVAT requires local install; myvariant is faster for batch annotation. Switch to OpenCRAVAT for PHI-sensitive or offline workflows. |

## References

- Lelong S et al. 2022. BioThings SDK: a toolkit for building high-performance data APIs in biomedical research. *Bioinformatics* 38:2077.
- Xin J et al. 2016. High-performance web services for querying gene and variant annotation. *Genome Biol* 17:91.
- Cheng J et al. 2023. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 381:eadg7492.
- Pejaver V et al. 2022. Calibration of computational tools for missense variant pathogenicity classification. *Am J Hum Genet* 109:2163.
- Pagel KA et al. 2020. Integrated informatics analysis of cancer-related variants. *JCO Clin Cancer Inform* 4:310. (OpenCRAVAT)
- Kopanos C et al. 2019. VarSome: the human genomic variant search engine. *Bioinformatics* 35:1978.
- Stawiński P, Płoski R. 2024. Genebe.net: implementation and validation of an automatic ACMG variant pathogenicity criteria assignment. *Clin Genet* 106:119.
- myvariant.info docs: `https://docs.myvariant.info/en/latest/`
- BioThings field metadata: `https://myvariant.info/v1/metadata/fields`

## Related Skills

- clinical-databases/clinvar-lookup - Source-level deep ClinVar queries
- clinical-databases/gnomad-frequencies - Source-level deep gnomAD queries
- clinical-databases/dbsnp-queries - Source-level rsID resolution
- clinical-databases/acmg-classification - ACMG framework with Pejaver calibration
- clinical-databases/variant-prioritization - Pipeline using aggregated annotations
