'''Fetch a ClinVar (VCV) or dbSNP record by UID and print a small dict. Both XML responses lack a usable DTD,
so Entrez.read() cannot parse them; this uses xml.etree.ElementTree. Report only the record's own stated
classification -- never a diagnosis or treatment recommendation.

Inputs: --email (required by NCBI), --db clinvar|snp, --uid (e.g. ClinVar UID 4887763, dbSNP UID 429358 for rs429358).
Usage: python variant_records.py --email you@inst.edu --db clinvar --uid 4887763
Import: from variant_records import clinvar_record, snp_record   (set Entrez.email first)
Tested: Biopython 1.88, live NCBI E-utilities 2026-09-21.
'''
import argparse
import xml.etree.ElementTree as ET
from Bio import Entrez


def clinvar_record(uid):
    h = Entrez.efetch(db='clinvar', id=uid, rettype='vcv', retmode='xml')
    raw = h.read(); h.close()
    text = raw.decode() if isinstance(raw, bytes) else raw
    archive = ET.fromstring(text).find('.//VariationArchive')
    sig = archive.find('.//Classifications/GermlineClassification/Description')
    return {
        'accession': archive.get('Accession'),
        'variation_name': archive.get('VariationName'),
        'clinical_significance': sig.text if sig is not None else None,
    }


def snp_record(uid):
    h = Entrez.efetch(db='snp', id=uid, rettype='xml', retmode='xml')
    raw = h.read(); h.close()
    text = raw.decode() if isinstance(raw, bytes) else raw
    ns = {'s': 'https://www.ncbi.nlm.nih.gov/SNP/docsum'}
    doc = ET.fromstring(text).find('s:DocumentSummary', ns)
    gene = doc.find('.//s:GENE_E/s:NAME', ns)
    sig = doc.find('s:CLINICAL_SIGNIFICANCE', ns)
    return {
        'chr': doc.findtext('s:CHR', default=None, namespaces=ns),
        'gene': gene.text if gene is not None else None,
        'clinical_significance': sig.text if sig is not None else None,
    }


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--email', required=True)
    ap.add_argument('--db', choices=['clinvar', 'snp'], required=True)
    ap.add_argument('--uid', required=True)
    a = ap.parse_args()
    Entrez.email = a.email
    print((clinvar_record if a.db == 'clinvar' else snp_record)(a.uid))
