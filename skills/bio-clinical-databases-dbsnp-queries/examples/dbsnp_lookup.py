'''dbSNP Build 156 query patterns: rsID resolution, merge-chain following, SPDI conversion, ALFA.

Reference: requests 2.31+, myvariant 1.0+ | checked against live NCBI Variation Services v0 on 2026-09-15.
Build 156 (Sept 2022) is the current schema; pre-Build-156 E-utilities returns thin legacy summary.
'''
import requests
import time
import myvariant
import pandas as pd

VARSVC = 'https://api.ncbi.nlm.nih.gov/variation/v0'

REFSEQ_GRCH38 = {
    '1': 'NC_000001.11', '2': 'NC_000002.12', '3': 'NC_000003.12', '4': 'NC_000004.12',
    '5': 'NC_000005.10', '6': 'NC_000006.12', '7': 'NC_000007.14', '8': 'NC_000008.11',
    '9': 'NC_000009.12', '10': 'NC_000010.11', '11': 'NC_000011.10', '12': 'NC_000012.12',
    '13': 'NC_000013.11', '14': 'NC_000014.9', '15': 'NC_000015.10', '16': 'NC_000016.10',
    '17': 'NC_000017.11', '18': 'NC_000018.10', '19': 'NC_000019.10', '20': 'NC_000020.11',
    '21': 'NC_000021.9', '22': 'NC_000022.11', 'X': 'NC_000023.11', 'Y': 'NC_000024.10'
}

ALFA_BIOPROJECT = 'PRJNA507278'
ALFA_BIOSAMPLES = {
    'SAMN10492705': 'Total', 'SAMN10492695': 'European', 'SAMN10492703': 'African',
    'SAMN10492696': 'African Others', 'SAMN10492698': 'African American', 'SAMN10492704': 'Asian',
    'SAMN10492697': 'East Asian', 'SAMN10492701': 'Other Asian', 'SAMN10492702': 'South Asian',
    'SAMN10492699': 'Latin American 1', 'SAMN10492700': 'Latin American 2', 'SAMN11605645': 'Other'
}


def refsnp(rsid, sleep=0.34):
    '''Fetch full Build 156 RefSNP JSON. sleep=0.34s -> ~3 req/s without API key.'''
    rs_int = str(rsid).lstrip('rs')
    r = requests.get(f'{VARSVC}/refsnp/{rs_int}', timeout=30)
    time.sleep(sleep)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()


def resolve_merge_chain(rsid, max_hops=10):
    '''Follow multi-hop merge chain. Cycle-safe with max_hops cap.

    A merged record has no primary_snapshot_data; merged_snapshot_data is an object
    {"merged_into": ["429358"], ...} (rs630496 -> rs429358). Each hop is one request.
    '''
    chain = []
    current = str(rsid).lstrip('rs')
    for _ in range(max_hops):
        if current in chain:
            return {'error': 'merge_cycle', 'chain': chain}
        chain.append(current)
        payload = refsnp(current)
        if payload is None:
            return {'status': 'not_found', 'final_rsid': current, 'chain': chain}
        if payload.get('is_withdrawn'):
            return {'status': 'withdrawn', 'final_rsid': current, 'chain': chain,
                    'reason': payload.get('withdrawn_release', {})}
        if payload.get('primary_snapshot_data') is not None:
            return {'status': 'resolved', 'final_rsid': payload.get('refsnp_id'), 'chain': chain,
                    'merged_from': [m.get('merged_rsid') for m in payload.get('dbsnp1_merges', [])]}
        merged_into = (payload.get('merged_snapshot_data') or {}).get('merged_into') or []
        if not merged_into:
            return {'status': 'orphan', 'final_rsid': current, 'chain': chain}
        current = str(merged_into[0])
    return {'error': 'hop_limit_exceeded', 'chain': chain}


def alleles_grch38(payload):
    '''Extract ALT alleles from RefSNP JSON for GRCh38 only. Returns list of {ref, alt, spdi}.

    More than one ALT allele = multi-allelic cluster (rs334: T>A, T>C, T>G).
    Assembly traits sit under placement_annot.seq_id_traits_by_assembly.
    '''
    if payload is None or payload.get('is_withdrawn'):
        return []
    primary = payload.get('primary_snapshot_data', {})
    out = []
    for placement in primary.get('placements_with_allele', []):
        traits = (placement.get('placement_annot') or {}).get('seq_id_traits_by_assembly') or [{}]
        if 'GRCh38' not in (traits[0].get('assembly_name') or ''):
            continue
        for allele in placement.get('alleles', []):
            spdi = allele.get('allele', {}).get('spdi', {})
            if spdi.get('inserted_sequence') == spdi.get('deleted_sequence'):
                continue
            out.append({
                'ref': spdi.get('deleted_sequence'),
                'alt': spdi.get('inserted_sequence'),
                'spdi': f"{spdi.get('seq_id')}:{spdi.get('position')}:{spdi.get('deleted_sequence')}:{spdi.get('inserted_sequence')}",
                'pos_0based': spdi.get('position')
            })
    return out


def vcf_to_canonical_spdi(chrom, pos, ref, alt):
    '''VCF (1-based, GRCh38) -> canonical SPDI string (0-based).

    The endpoint returns seq_id/position/deleted_sequence/inserted_sequence directly under data,
    plus data.warnings when REF does not match the reference (raised here).
    '''
    refseq = REFSEQ_GRCH38.get(str(chrom).removeprefix('chr'))
    if refseq is None:
        return None
    raw_spdi = f'{refseq}:{pos - 1}:{ref}:{alt}'
    r = requests.get(f'{VARSVC}/spdi/{raw_spdi}/canonical_representative', timeout=30)
    time.sleep(0.34)
    if not r.ok:
        return None
    d = r.json().get('data', {})
    if d.get('warnings'):
        raise ValueError(f"{chrom}:{pos} {ref}>{alt}: {[w.get('message') for w in d['warnings']]}")
    return f"{d['seq_id']}:{d['position']}:{d['deleted_sequence']}:{d['inserted_sequence']}"


def spdi_to_rsid(spdi_str):
    '''SPDI -> rsID if a cluster exists at that position with matching allele.'''
    r = requests.get(f'{VARSVC}/spdi/{spdi_str}/rsids', timeout=30)
    time.sleep(0.34)
    if not r.ok:
        return None
    rsids = r.json().get('data', {}).get('rsids', [])
    return rsids[0] if rsids else None


def hgvs_to_canonical_spdi(hgvs):
    '''HGVS-g / HGVS-c -> canonical SPDI via Variation Services contextuals.'''
    r = requests.get(f'{VARSVC}/hgvs/{hgvs}/contextuals', timeout=30)
    time.sleep(0.34)
    if not r.ok:
        return None
    contextuals = r.json().get('data', {}).get('spdis', [])
    return contextuals[0] if contextuals else None


def _clinvar_significance(entry):
    '''ClinVar significance from a myvariant record: clinvar.rcv (list or single dict).'''
    rcv = (entry.get('clinvar') or {}).get('rcv')
    rcvs = rcv if isinstance(rcv, list) else ([rcv] if rcv else [])
    sigs = sorted({x.get('clinical_significance') for x in rcvs if x.get('clinical_significance')})
    return ';'.join(sigs) or None


def batch_normalize_rsids(rsids):
    '''Normalize a list of rsIDs: resolve merges, flag multi-allelic, one row per input rsID.

    myvariant returns one hit per allele for a multi-allelic rsID (gnomAD values there are gnomAD 2.1.1);
    hits are collapsed per rsID with per-allele values joined by ';'.
    '''
    mv = myvariant.MyVariantInfo()
    fields = ['dbsnp.rsid', 'gnomad_exome.af.af', 'gnomad_genome.af.af', 'clinvar.rcv.clinical_significance']
    unique_rsids = list(dict.fromkeys(rsids))
    # Resolve merges first: myvariant.info's dbsnp collection is keyed by the CURRENT rsID, so a
    # merged input (e.g. rs630496) returns notfound unless queried by its resolved canonical rsID
    # (e.g. rs429358). Query by canonical rsID, not the input rsID. resolve_merge_chain's
    # final_rsid is a bare digit string (e.g. '429358'); re-add the 'rs' prefix for the query.
    merge_results = {rsid: resolve_merge_chain(rsid) for rsid in unique_rsids}
    def _rs(x):
        return x if str(x).startswith('rs') else f'rs{x}'
    canonical = {rsid: _rs(merge_results[rsid].get('final_rsid', rsid)) for rsid in unique_rsids}
    hits_by_rsid = {}
    for entry in mv.getvariants(list(dict.fromkeys(canonical.values())), fields=fields):
        hits_by_rsid.setdefault(entry.get('query'), []).append(entry)
    rows = []
    for rsid in unique_rsids:
        merge_result = merge_results[rsid]
        canonical_rsid = canonical[rsid]
        hits = [h for h in hits_by_rsid.get(canonical_rsid, []) if not h.get('notfound')]
        alleles = alleles_grch38(refsnp(canonical_rsid))
        # one ';'-separated value per myvariant hit, '.' where the hit has no value (keeps alleles aligned)
        join = lambda vals: ';'.join('.' if v is None else str(v) for v in vals) or None
        rows.append({
            'input_rsid': rsid,
            'canonical_rsid': merge_result.get('final_rsid'),
            'status': merge_result.get('status'),
            'chain_length': len(merge_result.get('chain', [])),
            'is_multiallelic': len(alleles) > 1,
            'grch38_alleles': ','.join(f"{a['ref']}>{a['alt']}" for a in alleles),
            'myvariant_ids': join(h.get('_id') for h in hits),
            'gnomad_v2_exome_af': join((h.get('gnomad_exome') or {}).get('af', {}).get('af') for h in hits),
            'gnomad_v2_genome_af': join((h.get('gnomad_genome') or {}).get('af', {}).get('af') for h in hits),
            'clinvar_sig': join(_clinvar_significance(h) for h in hits)
        })
    return pd.DataFrame(rows)


def alfa_population_frequencies(rsid):
    '''ALFA per-population allele counts from /refsnp/{id}/frequency (not in the RefSNP JSON).

    ALFA aggregates dbGaP studies across 12 ancestry groups. Use when:
    - Variant is array-genotyped (gnomAD may miss it)
    - Consent-respecting frequency lookup needed
    Do NOT use for rare-variant FAF95 -- use gnomAD instead.
    '''
    rs_int = str(rsid).lstrip('rs')
    r = requests.get(f'{VARSVC}/refsnp/{rs_int}/frequency', timeout=30)
    time.sleep(0.34)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    out = {}
    for placement, record in r.json().get('results', {}).items():
        counts = record.get('counts', {}).get(ALFA_BIOPROJECT, {}).get('allele_counts', {})
        for biosample, alleles in counts.items():
            total = sum(alleles.values())
            out[ALFA_BIOSAMPLES.get(biosample, biosample)] = {
                'placement': placement,
                'ref': record.get('ref'),
                'counts': alleles,
                'total': total,
                'af': {a: n / total for a, n in alleles.items()} if total else None
            }
    return out


if __name__ == '__main__':
    chain = resolve_merge_chain('rs630496')
    print(f'rs630496 (merged): status={chain["status"]}, final=rs{chain["final_rsid"]}, chain={chain["chain"]}')
    payload = refsnp(chain['final_rsid'])
    alleles = alleles_grch38(payload)
    print(f'GRCh38 alleles: {alleles}; merged_from: {chain.get("merged_from")}')

    spdi = vcf_to_canonical_spdi('17', 43106487, 'A', 'C')  # BRCA1 c.181T>G
    print(f'BRCA1 chr17:43106487:A>C canonical SPDI: {spdi}')
    print(f'  -> rsID: {spdi_to_rsid(spdi)}')

    alfa = alfa_population_frequencies('rs6025')
    print(f'rs6025 ALFA Total: {alfa.get("Total")}')
