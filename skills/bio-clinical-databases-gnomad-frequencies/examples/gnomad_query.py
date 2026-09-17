'''gnomAD v4 query patterns: GraphQL, grpmax FAF95, ACMG BS1/BA1 research tags, LOEUF constraint.

Reference: requests 2.31+, myvariant 1.0+, hail 0.2.130+ | checked against the live gnomAD GraphQL API on 2026-09-15.
v4.1 (May 2024) is current; v4.0 had AN under-counting fixed in v4.1.
Variants from patients or research participants sent to a public API can disclose them: use only under
the cohort's consent and institutional approvals.
'''
import time
import requests
import myvariant
import pandas as pd

GNOMAD_API = 'https://gnomad.broadinstitute.org/api'

BOTTLENECK_GROUPS = {'ami', 'asj', 'fin', 'remaining'}


def _post_graphql(query, variables, max_retries=5, base_delay=2.0):
    '''POST to the gnomAD GraphQL API with bounded exponential backoff on HTTP 429.

    The browser API rate-limits sustained querying: HTTP 429 (a plain HTML page, not JSON) was
    observed live after one gene-variant-list query plus a handful of single lookups (checked
    2026-09-15). For more than a few dozen variants, use the Hail Table snippet below instead of
    looping GraphQL calls.
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


def query_variant_v4(chrom, pos, ref, alt):
    '''Query a single GRCh38 variant in gnomAD v4 via GraphQL. Returns exome + genome data + grpmax FAF95.

    Coordinates MUST be GRCh38: GRCh37 ids return "Variant not found", the same answer as a truly
    absent variant. Returns None when not found; raises on any other GraphQL error.
    grpmax_faf95 is the ACMG-grade frequency: 95% lower-CI of AF, excluding bottleneck groups
    (AMI, ASJ, FIN, REMAINING) by design. Use this -- NOT raw AF -- for BS1/BA1 application.
    '''
    query = '''
    query VariantById($variantId: String!) {
      variant(variantId: $variantId, dataset: gnomad_r4) {
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
          filters
          populations { id ac an }
          faf95 { popmax popmax_population }
        }
      }
    }
    '''
    variant_id = f'{chrom}-{pos}-{ref}-{alt}'
    body = _post_graphql(query, {'variantId': variant_id})
    errors = [e.get('message') for e in body.get('errors') or []]
    if errors and errors != ['Variant not found']:
        raise RuntimeError(f'gnomAD GraphQL error for {variant_id}: {errors}')
    return (body.get('data') or {}).get('variant')


def grpmax_faf95(payload):
    '''Extract grpmax FAF95 -- the ACMG-grade frequency. Prefer exome over genome.

    faf95 is None for absent variants (source='absent') and for present variants whose FAF95 is
    undefined because too few alleles were seen (source='present_faf95_undefined').
    Bottleneck groups (AMI, ASJ, FIN, REMAINING) are excluded by design.
    '''
    if payload is None:
        return {'faf95': None, 'grpmax_ancestry': None, 'source': 'absent'}
    for source in ('exome', 'genome'):
        data = payload.get(source) or {}
        faf = data.get('faf95') or {}
        if faf.get('popmax') is not None:
            return {
                'faf95': faf['popmax'],
                'grpmax_ancestry': faf.get('popmax_population'),
                'source': source,
                'pass': 'PASS' in (data.get('filters') or []) or not data.get('filters')
            }
    return {'faf95': None, 'grpmax_ancestry': None, 'source': 'present_faf95_undefined'}


def max_credible_af(prevalence, max_allelic_contribution=1.0, max_genetic_contribution=1.0,
                    penetrance=1.0):
    '''Whiffin 2017 max-credible-AF formula for ACMG BS1 application.

    Args:
        prevalence: disease prevalence (e.g., HCM = 1/500 = 0.002)
        max_allelic_contribution: max share of cases attributable to a single allele
        max_genetic_contribution: max share of cases attributable to this gene
        penetrance: probability variant carriers develop disease

    Returns: gene-specific max-credible per-allele frequency under dominant inheritance.
             For autosomal recessive, multiply by carrier-frequency-squared appropriately.
    '''
    return (prevalence * max_genetic_contribution * max_allelic_contribution) / (penetrance * 2)


def apply_acmg_freq_codes(faf95_val, max_credible, ba1_threshold=0.05):
    '''Apply ClinGen SVI BS1/BA1/PM2_Supporting frequency tags from grpmax_faf95 (research annotation).

    BA1 (stand-alone benign): default 5% in non-bottleneck group. VCEPs override
    (Hearing Loss = 0.5% AR). Check the relevant VCEP CSpec.
    BS1: gene-specific max-credible-AF via Whiffin formula.
    PM2_Supporting: absent or ultra-rare (faf95_val None = absent or FAF95 undefined).
    '''
    if faf95_val is None:
        return 'PM2_Supporting'
    if faf95_val > ba1_threshold:
        return 'BA1'
    if faf95_val > max_credible:
        return 'BS1'
    return None


def query_gene_constraint_v4(gene_symbol):
    '''Pull gene constraint from v4 (chrX genes included). Falls back only when v4 returns none.'''
    query = '''
    query GeneConstraint($symbol: String!) {
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
        gene['note'] = 'no v4 constraint returned -- try v2.1.1 (reference_genome: GRCh37)'
    loeuf = (gene.get('gnomad_constraint') or {}).get('oe_lof_upper')
    gene['loeuf_v4_first_decile'] = loeuf is not None and loeuf < 0.6
    return gene


def annotate_variant_list(variants):
    '''Batch gnomAD AF via the myvariant.info aggregator.

    myvariant.info carries gnomAD 2.1.1 (GRCh37) and no FAF95 fields (metadata, 2026-09-15), so this
    returns raw v2.1.1 exome/genome AF only. Use query_variant_v4() for v4 grpmax FAF95.
    '''
    mv = myvariant.MyVariantInfo()
    fields = ['gnomad_exome.af.af', 'gnomad_exome.an.an', 'gnomad_genome.af.af']
    results = mv.getvariants(variants, fields=fields)
    rows = []
    for r in results:
        exome = r.get('gnomad_exome') or {}
        genome = r.get('gnomad_genome') or {}
        rows.append({
            'variant': r.get('query'),
            'myvariant_id': r.get('_id'),
            'exome_af_v2': (exome.get('af') or {}).get('af'),
            'exome_an_v2': (exome.get('an') or {}).get('an'),
            'genome_af_v2': (genome.get('af') or {}).get('af')
        })
    return pd.DataFrame(rows)


def hail_bulk_filter_snippet():
    '''Print Hail code for bulk rare-variant filtering against gnomAD v4 exomes.

    Reference: hail 0.2.130+. Requires GCS authentication.
    Use when filtering >10k variants -- API rate limits make GraphQL impractical at scale.
    Patient exomes: run only in an approved environment under the cohort's approvals.
    '''
    return '''
import hail as hl
hl.init(default_reference='GRCh38')

ht_v4 = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.sites.ht')
mt = hl.import_vcf('input.vcf.gz', force_bgz=True, reference_genome='GRCh38')
mt = mt.annotate_rows(gnomad=ht_v4[mt.locus, mt.alleles])

# Filter to variants below grpmax FAF95 threshold; absent variants pass
mt = mt.filter_rows(
    (hl.is_missing(mt.gnomad)) |
    (hl.is_missing(mt.gnomad.grpmax_faf95)) |
    (mt.gnomad.grpmax_faf95.faf95 < 0.0001)
)
mt.write('rare_variants.mt', overwrite=True)
'''


if __name__ == '__main__':
    brca1_variant = query_variant_v4('17', 43106487, 'A', 'C')  # BRCA1 c.181T>G, GRCh38
    grpmax_info = grpmax_faf95(brca1_variant)
    print(f'BRCA1 c.181T>G grpmax FAF95: {grpmax_info}')

    hcm_threshold = max_credible_af(prevalence=1/500, max_genetic_contribution=0.30,
                                     max_allelic_contribution=0.10, penetrance=0.80)
    print(f'HCM gene-specific max-credible-AF (MYH7-like): {hcm_threshold:.6f}')
    code = apply_acmg_freq_codes(grpmax_info['faf95'], hcm_threshold)
    print(f'Frequency tag (research annotation): {code}')

    for symbol in ('SCN2A', 'DMD'):
        constraint = query_gene_constraint_v4(symbol)
        if constraint:
            print(f'{symbol} (chr{constraint["chrom"]}) LOEUF (v4): '
                  f'{(constraint.get("gnomad_constraint") or {}).get("oe_lof_upper")}; '
                  f'first decile: {constraint["loeuf_v4_first_decile"]}')
