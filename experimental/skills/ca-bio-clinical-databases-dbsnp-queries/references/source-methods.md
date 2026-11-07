# Historical source methods

Source build, service and frequency claims have not been reverified; check the actual current contract and release. Inspect all allele mappings rather than assuming uniqueness.

## rsID Is a Cluster Identifier, Not a Variant Identifier

This is the load-bearing concept. dbSNP cluster definition: ss records (submitted SNPs) are mapped to the genome and clustered into RefSNPs by *position + variant type*, not by allele. A single rsID can point to a *locus* with multiple alleles:

- `rs12345` may resolve to {A>G, A>T, A>C} at one position; the RefSNP JSON `primary_snapshot_data.placements_with_allele[*].alleles` enumerates them.
- ~6-8% of dbSNP rsIDs are multi-allelic.
- PLINK and many older tools historically misuse rsIDs as if they were variant identifiers, which fails for multi-allelic sites and yields wrong genotype assignments.

**Rule:** Use rsID as a *human-facing label only*; use SPDI or ClinGen Allele Registry CA ID for joins.

## Build 156 Schema Overhaul: What Changed

| Aspect | Build 151 (2017) | Build 156 (2022) and current |
|--------|------------------|------------------------------|
| Distribution | Relational SQL dumps + XML | JSON per RefSNP, partitioned by chromosome |
| FTP path | `ftp/snp/organisms/human_9606/` | `ftp.ncbi.nlm.nih.gov/snp/latest_release/JSON/` |
| Primary key | `snp_id`, `ss_id` | `refsnp_id`, with `primary_snapshot_data` block |
| Frequency data | Embedded sparse | ALFA aggregated populations |
| Merge tracking | `RsMergeArch.bcp.gz` | `refsnp-merged.json.bz2` (also `RsMergeArch.bcp.gz` retained for legacy) |
| Withdrawn | `SNPHistory.bcp.gz` | `refsnp-withdrawn.json.bz2` |
| API access | Legacy E-utilities `db=snp` only | Variation Services REST `/v0/refsnp/{id}` returns the full JSON |

E-utilities still works for `db=snp` but returns a thin pre-156 summary missing key fields like `primary_snapshot_data.placements_with_allele`; pipelines reliant on Entrez get out-of-date data.

## RsMergeArch: The Multi-Hop Merge Footgun

When two rsIDs are found to refer to the same allele cluster, the higher (later-assigned) rsID is merged into the lower. `RsMergeArch.bcp.gz` stores `(rsHigh, rsLow, rsCurrent)` tuples.

**The trap:** `rsCurrent` in any given row is the merge target at the time of that merge event, NOT the current dbSNP rsID. A multi-merge chain (rs3 -> rs2 -> rs1, then later rs1 -> rs0) appears as multiple rows. Naive one-hop lookup resolves to a stale ID.

Withdrawn rsIDs (submitter-withdrawn or QC-failed) live in `SNPHistory.bcp.gz`, not RsMergeArch. Both tables must be consulted to resolve any historical rsID.

## SPDI: The Canonical Variant Representation

SPDI (Sequence:Position:Deletion:Insertion) format: `NC_000017.11:43044294:G:A`. Position is **0-based, half-open** (differs from HGVS's 1-based, fully-closed).

The **Contextual Allele** transformation (Variant Overprecision Correction Algorithm) returns the right-aligned, normalized canonical form across left-aligned VCFs and right-aligned HGVS conventions. This is the basis for ClinGen Allele Registry CA ID computation.

| Representation | Build dependency | Transcript dependency | Bijective? | Best for |
|----------------|-----------------|----------------------|------------|----------|
| **VCF (chrom-pos-ref-alt)** | Yes | No | Yes (same build) | Pipelines, bulk |
| **SPDI** | Yes (via RefSeq accession) | No | Yes for SNV/small indel | Canonical normalization |
| **HGVS-g** | Yes (`NC_xxxxx.N`) | No | Yes for SNV/small indel | Human-readable genomic |
| **HGVS-c** | Indirect (via transcript) | Yes | No (one HGVS-c -> many HGVS-g) | Clinical reporting |
| **HGVS-p** | Indirect | Yes | Degenerate (one HGVS-p -> many HGVS-c) | Protein-level annotation |
| **rsID** | None (cluster identifier) | None | NO (multi-allelic) | Human label only |
| **CA ID** | None (canonical) | None | Yes | Cross-database join |

## Decision Tree by Query Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Resolve single rsID to coordinates + alleles | Variation Services `/v0/refsnp/{id}` | Returns full Build 156 JSON, including merge history |
| Resolve historical/deprecated rsID | Variation Services `/v0/refsnp/{id}` -> follow `merged_snapshot_data` chain | Single-hop RsMergeArch lookup misses multi-hop chains |
| Batch query 100-10k rsIDs | myvariant.info `getvariants(rsids)` | Aggregated with ClinVar/gnomAD overlay; rate-limit safe |
| Convert coords <-> rsID | myvariant.info HGVS query or Variation Services `/spdi/{spdi}/rsid` | SPDI is the canonical bridge |
| Normalize variant representations | Variation Services `/hgvs/{hgvs}/contextuals` | Returns canonical SPDI, right-aligned |
| Bulk genomic-wide rsID -> coords | Local download of `refsnp-chr{N}.json.bz2` + parser | No rate limits; weekly snapshots |
| Joining dbSNP with gnomAD by ID | Use SPDI or CA ID, never rsID alone | rsID is a cluster; alleles may not match |
| Get population AF for common variant | ALFA (via Variation Services) for array-genotyped variants; gnomAD for sequencing-derived | Different sample compositions |

## ALFA Frequencies vs gnomAD

| Source | Sample basis | Variants covered | When to use |
|--------|-------------|------------------|-------------|
| **ALFA** | ~1M dbGaP subjects (array + WGS) across 12 ancestry groups | 447M+ sites; broader (includes array-only) | Common variants, dbGaP-deposited cohorts |
| **gnomAD v4** | 807k WGS+WES individuals across 9 ancestry groups | Sequencing-derived (deeper at rare variants) | Rare-variant FAF95, ACMG BS1/BA1 |

ALFA does NOT provide FAF95-style upper-bound CIs; raw AF only. ALFA captures consent-tier metadata enabling consent-respecting lookups for variants gnomAD doesn't carry.

ALFA counts are NOT in the RefSNP JSON (its `frequency` list holds 1000 Genomes, gnomAD, TOPMed and other studies). They come from `GET /variation/v0/refsnp/{id}/frequency`, under BioProject `PRJNA507278`, keyed by ALFA BioSample ID.

```python
ALFA_BIOPROJECT = 'PRJNA507278'
ALFA_BIOSAMPLES = {
    'SAMN10492705': 'Total', 'SAMN10492695': 'European', 'SAMN10492703': 'African',
    'SAMN10492696': 'African Others', 'SAMN10492698': 'African American', 'SAMN10492704': 'Asian',
    'SAMN10492697': 'East Asian', 'SAMN10492701': 'Other Asian', 'SAMN10492702': 'South Asian',
    'SAMN10492699': 'Latin American 1', 'SAMN10492700': 'Latin American 2', 'SAMN11605645': 'Other'
}

def alfa_frequency(rsid, population='Total'):
    '''ALFA allele counts and frequencies for one population from /refsnp/{id}/frequency.

    Returns {placement: {'ref', 'total_alleles', 'af': {allele: freq}}}; placement keys look like
    '1@169549810' for rs6025.
    '''
    rs_int = str(rsid).lstrip('rs')
    r = requests.get(f'{VARSVC}/refsnp/{rs_int}/frequency', timeout=30)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    out = {}
    for placement, record in r.json().get('results', {}).items():
        counts = record.get('counts', {}).get(ALFA_BIOPROJECT, {}).get('allele_counts', {})
        for biosample, alleles in counts.items():
            if ALFA_BIOSAMPLES.get(biosample) != population:
                continue
            total = sum(alleles.values())
            out[placement] = {'ref': record.get('ref'), 'total_alleles': total,
                              'af': {a: n / total for a, n in alleles.items()} if total else None}
    return out or None
```

## Per-Operation Failure Modes

**1. Treating rsID as a unique variant identifier**
- Trigger: Join two databases by rsID expecting a single variant.
- Mechanism: rsID is a cluster identifier; multi-allelic clusters have 2-4 alleles at one position.
- Symptom: Allele mismatches at low rate (~6-8% of sites); silent merger of unrelated variants.
- Fix: Normalize both sides to SPDI or CA ID before joining.

**2. Single-hop RsMergeArch lookup**
- Trigger: Read one row of `RsMergeArch.bcp.gz` and treat `rsCurrent` as the final answer.
- Mechanism: Multi-hop merges (rs3 -> rs2 -> rs1 -> rs0) span multiple rows; each row records one hop only.
- Symptom: Resolved rsID is itself stale; subsequent queries return outdated annotation.
- Fix: Follow merge chains via Variation Services: each merged record's `merged_snapshot_data.merged_into`, one hop per request, until a record with `primary_snapshot_data`. Current records list the ids they absorbed in `dbsnp1_merges`.

**3. Confusing withdrawn vs merged**
- Trigger: Query a withdrawn rsID and find no merge target.
- Mechanism: Withdrawn rsIDs live in `SNPHistory.bcp.gz`, not RsMergeArch. They are NOT merged into anything; the cluster was QC-failed or submitter-retracted.
- Symptom: 404 from naive lookup; pipelines proceed with stale rsID.
- Fix: Check `is_withdrawn` field in RefSNP JSON; flag the variant for manual review.

**4. E-utilities thin summary**
- Trigger: Use `Entrez.esummary(db='snp')` and treat output as authoritative.
- Mechanism: E-utilities `db=snp` returns a pre-Build-156 summary missing `primary_snapshot_data.placements_with_allele`, frequency data, and merge history.
- Symptom: Annotations look incomplete or out-of-date.
- Fix: Use Variation Services REST `/v0/refsnp/{id}` for full JSON.

**5. SPDI 0-based vs HGVS 1-based mismatch**
- Trigger: Build SPDI string from VCF position without converting to 0-based.
- Mechanism: SPDI position is 0-based half-open; VCF is 1-based.
- Symptom: Coordinate off-by-one; SPDI does not resolve to expected rsID.
- Fix: SPDI position = (VCF position - 1). For indels, also normalize ref/alt.

**6. Strand/orientation ambiguity for A/T C/G variants**
- Trigger: Merge variants across builds or platforms relying on rsID alone.
- Mechanism: rsID is locus-level; opposite-strand alleles get the same rsID with different ref/alt representation.
- Symptom: Strand-flipped genotypes after merge.
- Fix: Use SPDI (which encodes strand via deleted/inserted sequence) or MAF-match for ambiguous variants.

## Reconciliation: When Sources Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| dbSNP rsID returns 404 in current build | Withdrawn (in `refsnp-withdrawn.json.bz2`) | Check withdrawal reason; consider manual curation |
| rsID resolves to different coords across builds | Genome assembly change (GRCh37 -> GRCh38) | Use SPDI with explicit RefSeq accession; lift over via `pyliftover` |
| ALFA AF and gnomAD AF disagree by >2x | Different sample compositions; ALFA includes array-only sites under-represented in gnomAD | Trust gnomAD for sequencing data; ALFA for array-derived; use the more relevant source |
| Multiple rsIDs map to one SPDI | True duplicates from independent submissions; rare since Build 152 enforced cluster merging | Pick the lowest rsID per RsMergeArch convention |
| One rsID has different ref allele in dbSNP vs gnomAD | dbSNP uses NCBI's ref; gnomAD aligns to its build assembly | Normalize to SPDI before joining |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| Build 156 | Current as of Sep 2022; JSON-only distribution | dbSNP NCBI |
| Multi-allelic rate | ~6-8% of dbSNP rsIDs are multi-allelic | operational estimate |
| Variation Services rate limit | 10 req/s with API key; 3 req/s without | NCBI E-utilities policy |
| SPDI position | 0-based, half-open | NCBI SPDI specification |
| HGVS position | 1-based, fully-closed | HGVS nomenclature |
| ALFA samples | ~1M individuals across 12 ancestries (2024 release) | NCBI dbGaP aggregation |
| Bulk download chunks | One file per chromosome (`refsnp-chr{N}.json.bz2`) | NCBI FTP |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| 404 from Variation Services on valid rsID | rsID is withdrawn or never assigned | Check `refsnp-withdrawn.json.bz2`; consider strand-flipped equivalent |
| Merge chain resolves but final rsID has different alleles | Multi-allelic cluster; pick allele matching the variant | Filter `placements_with_allele.alleles[*]` by allele match |
| ALFA frequency missing for common variant | Variant not in dbGaP-deposited studies | Fall back to gnomAD (sequencing-derived) |
| `Entrez.esummary` returns old data | Legacy E-utilities, not Build 156 schema | Switch to Variation Services REST `/v0/refsnp/{id}` |
| HGVS-c conversion fails for synonymous variants | Some HGVS-c rely on non-MANE transcripts not in NCBI default | Specify transcript explicitly; use VEP `--mane_select` |
| SPDI for indel does not round-trip | Left/right alignment mismatch | Use `/spdi/{spdi}/canonical_representative` for normalization |
| Bulk JSON parse OOM | `refsnp-chr1.json.bz2` is ~20GB uncompressed | Stream parse with `bz2.BZ2File` + line-by-line JSON; do not load whole file |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why not just use rsID for the join?" | rsID is a cluster identifier; ~6-8% of clusters are multi-allelic, causing silent mismatches. We use SPDI / CA ID. |
| "This annotation says rs12345 but the literature says rs67890" | rsIDs are merged; we resolved through `RsMergeArch` / `merged_snapshot_data` to the current canonical rsID. |
| "dbSNP frequency != gnomAD frequency" | ALFA (dbSNP-embedded) and gnomAD use different sample sets and different ascertainment (array vs sequencing); reconciled per use case. |
| "Why wasn't Entrez used?" | Entrez `db=snp` returns the pre-Build-156 summary missing key fields; we use Variation Services REST for the full JSON. |
| "Coordinate off-by-one in the SPDI" | SPDI is 0-based half-open; VCF is 1-based; intentional conversion applied. |

## References

- Phan L et al. 2025. The evolution of dbSNP: 25 years of impact in genomic research. *Nucleic Acids Res* 53:D925.
- Sayers EW et al. 2024. Database resources of the National Center for Biotechnology Information. *Nucleic Acids Res* 52:D33.
- Holmes JB et al. 2020. SPDI: data model for variants and applications at NCBI. *Bioinformatics* 36:1902.
- NCBI Variation Services API: `https://api.ncbi.nlm.nih.gov/variation/v0/`
- dbSNP FTP layout: `ftp.ncbi.nlm.nih.gov/snp/latest_release/JSON/`
- ALFA release notes: `https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/`
- ClinGen Allele Registry: `https://reg.clinicalgenome.org/docs/cg-car/`
