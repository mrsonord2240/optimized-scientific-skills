'''ClinVar query patterns: REST API, local VCF, ClinGen Allele Registry.

Reference: requests 2.31+, cyvcf2 0.30+ | checked against live E-utilities and the Allele Registry on 2026-09-15.
Handles 2024 v2 XML schema (VariationArchive anchor; tripartite Germline/Somatic/Oncogenicity).
ClinVar assertions are research evidence, not a diagnosis for an individual.
'''
import requests
import time
from cyvcf2 import VCF
import pandas as pd

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
ALLELE_REG = 'https://reg.clinicalgenome.org'

REVIEW_STATUS_STARS = {
    'practice guideline': 4,
    'reviewed by expert panel': 3,
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, single submitter': 1,
    'criteria provided, conflicting classifications': 1,
    'criteria provided, conflicting interpretations': 1,  # pre-2024 wording
    'no assertion criteria provided': 0,
    'no classification provided': 0
}


def clinvar_summary(variation_id):
    '''Fetch VCV-level summary by ClinVar VariationID via E-utilities esummary.'''
    r = requests.get(f'{EUTILS}/esummary.fcgi',
                     params={'db': 'clinvar', 'id': variation_id, 'retmode': 'json'},
                     timeout=30)
    r.raise_for_status()
    payload = r.json()['result'].get(str(variation_id), {})
    germline = payload.get('germline_classification', {})
    return {
        'vcv': payload.get('accession'),
        'name': payload.get('title'),
        'germline_class': germline.get('description'),
        'germline_review_status': germline.get('review_status'),
        'germline_last_evaluated': germline.get('last_evaluated'),
        'somatic_clinical_impact': payload.get('clinical_impact_classification', {}).get('description'),
        'oncogenicity': payload.get('oncogenicity_classification', {}).get('description'),
        'star_rating': REVIEW_STATUS_STARS.get(germline.get('review_status', '').lower(), 0)
    }


def car_record(hgvs_g):
    '''Fetch the ClinGen Allele Registry record for an HGVS-g with the public GET.

    PUT /allele needs authorization (HTTP 403). A REF that does not match the reference returns
    HTTP 400 IncorrectReferenceAllele; it is raised so it is not mistaken for "not registered".
    '''
    r = requests.get(f'{ALLELE_REG}/allele', params={'hgvs': hgvs_g}, timeout=30)
    if r.status_code == 400:
        err = r.json()
        raise ValueError(f"{err.get('errorType')}: {err.get('message')} (actual allele {err.get('actualAllele')})")
    r.raise_for_status()
    return r.json()


def car_id(hgvs_g):
    '''Resolve HGVS-g to ClinGen Allele Registry CA ID. Build/transcript-agnostic.'''
    at_id = car_record(hgvs_g).get('@id', '')
    return at_id.rsplit('/', 1)[-1] if at_id else None


def lookup_local_vcf(clinvar_vcf, chrom, pos, ref, alt):
    '''Query a local clinvar.vcf.gz by GRCh38 coords. Returns VCV-level INFO -- NOT condition-specific.

    ClinVar contigs are 1..22, X, Y, MT (no "chr"). For condition-stratified analysis use RCV-level XML parsing.
    '''
    vcf = VCF(clinvar_vcf)
    for v in vcf(f'{chrom}:{pos}-{pos}'):
        if v.REF == ref and alt in v.ALT:
            info = v.INFO
            return {
                'variation_id': v.ID,
                'allele_id': info.get('ALLELEID'),
                'clnsig': info.get('CLNSIG'),
                'clnsig_conf': info.get('CLNSIGCONF'),
                'clnrevstat': info.get('CLNREVSTAT'),
                'clndn': info.get('CLNDN'),
                'clnvc': info.get('CLNVC'),
                'clnhgvs': info.get('CLNHGVS'),
                'oncdn': info.get('ONCDN'),
                'scidn': info.get('SCIDN')
            }
    return None


def batch_resolve_to_car_then_clinvar(hgvs_list, sleep=0.34):
    '''Resolve HGVS-g list to CA IDs, then query ClinVar via the linked VariationID.

    sleep=0.34s -> ~3 req/s to stay under default NCBI rate limit without API key.
    '''
    rows = []
    for hgvs in hgvs_list:
        try:
            record = car_record(hgvs)
        except ValueError as err:
            rows.append({'hgvs': hgvs, 'ca_id': None, 'variation_id': None, 'error': str(err)})
            time.sleep(sleep)
            continue
        except requests.HTTPError as err:
            # Non-400 Registry errors (e.g. HTTP 500 "Unknown reference: <accession>" for an
            # unrecognized RefSeq accession) reach here; car_record only special-cases 400.
            # Record the Registry's own message instead of aborting the whole batch.
            try:
                message = err.response.json().get('message', str(err))
            except ValueError:
                message = str(err)
            rows.append({'hgvs': hgvs, 'ca_id': None, 'variation_id': None, 'error': message})
            time.sleep(sleep)
            continue
        time.sleep(sleep)
        at_id = record.get('@id', '')
        ca = at_id.rsplit('/', 1)[-1] if at_id else None
        clinvar_records = record.get('externalRecords', {}).get('ClinVarVariations', [])
        variation_id = clinvar_records[0]['variationId'] if clinvar_records else None
        row = {'hgvs': hgvs, 'ca_id': ca, 'variation_id': variation_id}
        if variation_id is not None:
            row.update(clinvar_summary(variation_id))
            time.sleep(sleep)
        rows.append(row)
    return pd.DataFrame(rows)


def filter_by_star_and_freshness(df, min_star=2, max_age_months=36):
    '''Keep variants with star >= min_star AND last_evaluated within max_age_months (research curation filter).

    Star=3 (VCEP) supersedes all lower-star records. The age rule also drops old expert-panel
    curations; exempt star 3 explicitly if that is not intended.
    Freshness rationale: Yauy 2022 documents ~1247 classification changes per monthly release.
    '''
    df = df.copy()
    df['last_evaluated_date'] = pd.to_datetime(df['germline_last_evaluated'], errors='coerce')
    today = pd.Timestamp.today()
    df['age_months'] = (today - df['last_evaluated_date']).dt.days / 30.44
    keep = (df['star_rating'] >= min_star) & (df['age_months'] <= max_age_months)
    return df[keep]


def parse_clnsig_conflict(clnsig_conf):
    '''Split CLNSIGCONF into individual conflicting calls.

    The meaning of "Conflicting classifications" depends on which calls conflict:
    - P vs LP : usually immaterial
    - P vs VUS : meaningful
    - P vs B/LB : major conflict; requires resolution
    '''
    if clnsig_conf is None:
        return {'calls': [], 'severity': None}
    calls = []
    for entry in str(clnsig_conf).split('|'):
        if '(' in entry:
            label, count = entry.rsplit('(', 1)
            calls.append({'label': label.strip('_'), 'submitter_count': int(count.rstrip(')'))})
    pathogenic_codes = {'Pathogenic', 'Likely_pathogenic'}
    benign_codes = {'Benign', 'Likely_benign'}
    has_path = any(c['label'] in pathogenic_codes for c in calls)
    has_benign = any(c['label'] in benign_codes for c in calls)
    has_vus = any(c['label'] == 'Uncertain_significance' for c in calls)
    if has_path and has_benign:
        conflict_severity = 'severe'
    elif has_path and has_vus:
        conflict_severity = 'meaningful'
    elif has_path:
        conflict_severity = 'minor'
    else:
        conflict_severity = 'non_pathogenic_only'
    return {'calls': calls, 'severity': conflict_severity}


def annotate_vcf_with_bcftools(input_vcf, clinvar_vcf, output_vcf):
    '''Emit the bcftools command to annotate a user VCF with 2024-schema INFO fields.

    input_vcf must use ClinVar contig names (17, not chr17), or every record annotates to ".".
    '''
    return (
        f'bcftools annotate '
        f'-a {clinvar_vcf} '
        f'-c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVC,INFO/CLNHGVS,'
        f'INFO/CLNSIGCONF,INFO/ONCDN,INFO/SCIDN '
        f'{input_vcf} -O z -o {output_vcf}'
    )


if __name__ == '__main__':
    example_hgvs = 'NC_000017.11:g.43106487A>C'  # BRCA1 c.181T>G, VCV000017661
    ca = car_id(example_hgvs)
    print(f'CA ID for {example_hgvs}: {ca}')
    summary = clinvar_summary(17661)
    print(f'Germline classification: {summary["germline_class"]} ({summary["star_rating"]} stars)')
    print(f'Last evaluated: {summary["germline_last_evaluated"]}')
