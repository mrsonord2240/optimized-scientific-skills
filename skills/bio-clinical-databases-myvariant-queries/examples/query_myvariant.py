'''myvariant.info aggregated annotation: batch queries, Lucene search, reproducibility.

Reference: myvariant 1.0+, requests 2.31+ | field paths and versions checked against the live service on 2026-09-15.
Records carry no _meta: per-source versions come from /v1/metadata (dbNSFP version is the
dominant staleness vector for in-silico predictors). myvariant carries gnomAD 2.1.1 (no FAF95)
and ClinVar 2025-05; use the gnomad-frequencies / clinvar-lookup Skills for current releases.
'''
import myvariant
import pandas as pd
import requests

mv = myvariant.MyVariantInfo()

CLINICAL_FIELDS = [
    'clinvar.rcv.clinical_significance',
    'clinvar.rcv.review_status',
    'clinvar.variant_id',
    'gnomad_exome.af.af',
    'gnomad_exome.an.an',
    'gnomad_genome.af.af',
    'dbsnp.rsid',
    'dbnsfp.alphamissense.score',
    'dbnsfp.alphamissense.pred',
    'dbnsfp.revel.score',
    'cadd.phred',
    'dbnsfp.bayesdel.add_af.score',
    'cosmic.cosmic_id',
    'civic.openCravatUrl'
]


def _as_list(value):
    return value if isinstance(value, list) else ([value] if value is not None else [])


def annotate_variants(hgvs_list):
    '''Batch-annotate a list of variants. Returns (DataFrame, versions_dict).

    getvariants() chunks to <=1000 ids per POST itself (a raw POST over 1000 returns HTTP 400).
    Versions come from /v1/metadata, once per run.
    '''
    rows = [_parse_record(r) for r in mv.getvariants(hgvs_list, fields=CLINICAL_FIELDS)]
    return pd.DataFrame(rows), metadata_versions()


def _parse_record(r):
    '''Defensively parse a myvariant record into a flat dict.'''
    if r.get('notfound'):
        return {'variant': r.get('query'), 'notfound': True}
    rcv = _as_list((r.get('clinvar') or {}).get('rcv'))
    gnomad_e = r.get('gnomad_exome') or {}
    gnomad_g = r.get('gnomad_genome') or {}
    dbnsfp = r.get('dbnsfp') or {}
    am = dbnsfp.get('alphamissense') or {}
    revel = dbnsfp.get('revel') or {}
    bayesdel = (dbnsfp.get('bayesdel') or {}).get('add_af') or {}
    return {
        'variant': r.get('query'),
        'myvariant_id': r.get('_id'),
        'rsid': (r.get('dbsnp') or {}).get('rsid'),
        'clinvar_sig': sorted({x.get('clinical_significance') for x in rcv if x.get('clinical_significance')}),
        'clinvar_review': sorted({x.get('review_status') for x in rcv if x.get('review_status')}),
        'clinvar_var_id': (r.get('clinvar') or {}).get('variant_id'),
        'gnomad_v2_af_exome': (gnomad_e.get('af') or {}).get('af'),
        'gnomad_v2_af_genome': (gnomad_g.get('af') or {}).get('af'),
        'alphamissense': am.get('score'),
        'alphamissense_pred': am.get('pred'),
        'revel': revel.get('score'),
        'cadd_phred': (r.get('cadd') or {}).get('phred'),
        'bayesdel_score': bayesdel.get('score'),
        'cosmic_id': (r.get('cosmic') or {}).get('cosmic_id')
    }


def find_pathogenic_in_gene(gene_symbol, max_results=500, min_star=2):
    '''Lucene search: ClinVar P/LP in a gene with star >= min_star. Returns the raw result dict;
    compare len(result['hits']) with result['total'] (hits holds at most max_results).

    Star convention (clinvar.rcv.review_status):
      "reviewed by expert panel" -> 3 stars (VCEP)
      "criteria provided, multiple submitters, no conflicts" -> 2 stars
      "criteria provided, single submitter" -> 1 star
    A record matches if any of its RCVs matches each term.
    '''
    review_query = ''
    if min_star >= 3:
        review_query = ' AND clinvar.rcv.review_status:"reviewed by expert panel"'
    elif min_star >= 2:
        review_query = (' AND (clinvar.rcv.review_status:"reviewed by expert panel" '
                        'OR clinvar.rcv.review_status:"criteria provided, multiple submitters, no conflicts")')
    q = (f'clinvar.gene.symbol:{gene_symbol} AND '
         'clinvar.rcv.clinical_significance:(Pathogenic OR "Likely pathogenic")' + review_query)
    return mv.query(q, size=max_results, fields=['_id', 'clinvar.rcv.clinical_significance',
                                                   'clinvar.rcv.review_status', 'clinvar.hgvs'])


def find_high_impact_in_region(chrom, start, end, min_cadd=25, am_threshold=0.564):
    '''Find variants in a GRCh37 (hg19) region with high CADD or high AlphaMissense.

    start/end are hg19 positions (hg38.start with cadd.phred returned 0 hits on 2026-09-15).
    Note: AlphaMissense developer threshold 0.564 is NOT the Pejaver 2022 calibrated
    PP3 threshold. Treat as supporting evidence only until ClinGen calibrates.
    '''
    q = (f'chrom:{chrom} AND hg19.start:[{start} TO {end}] AND '
         f'(cadd.phred:>{min_cadd} OR dbnsfp.alphamissense.score:>{am_threshold})')
    return mv.query(q, size=500, fields=['_id', 'cadd.phred', 'dbnsfp.alphamissense',
                                          'clinvar.rcv.clinical_significance'])


def metadata_versions():
    '''Pull current per-source versions across the entire myvariant.info instance.'''
    r = requests.get('https://myvariant.info/v1/metadata', timeout=30)
    r.raise_for_status()
    meta = r.json()
    return {name: src.get('version', src.get('build', '?'))
            for name, src in meta.get('src', {}).items()}


if __name__ == '__main__':
    print('Current myvariant.info per-source versions:')
    versions = metadata_versions()
    for src in ('clinvar', 'dbnsfp', 'gnomad', 'dbsnp', 'cosmic', 'civic'):
        print(f'  {src}: {versions.get(src)}')

    print('\nBatch annotating 4 test variants...')
    test = ['chr7:g.140453136A>T', 'rs121913529', 'rs1800566', 'rs104894155']
    df, src_versions = annotate_variants(test)
    print(df[['variant', 'rsid', 'clinvar_sig', 'gnomad_v2_af_exome',
              'alphamissense', 'revel', 'cadd_phred']].to_string(index=False))
    print(f'\ndbNSFP version recorded: {(src_versions or {}).get("dbnsfp", "n/a")}')
