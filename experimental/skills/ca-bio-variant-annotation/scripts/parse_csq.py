#!/usr/bin/env python
'''Flatten VEP CSQ or SnpEff ANN transcript blocks of an annotated VCF into one TSV row per block.

Inputs:  annotated VCF (plain or bgzip) carrying a VEP CSQ or SnpEff ANN INFO field.
Usage:   python scripts/parse_csq.py annotated.vcf.gz [--mane-only] [--impact HIGH,MODERATE] [--out parsed_annotations.tsv]
Output:  TSV with chrom, pos, ref, alt, then every field of the CSQ format line (or the 15 standard ANN fields).
         --mane-only keeps blocks with a populated MANE_SELECT (VEP, run with --mane_select); --impact keeps the
         listed IMPACT values. Impact is a triage bucket, not ACMG evidence.
Refuses: no CSQ or ANN header; --mane-only on a VCF whose CSQ format lacks MANE_SELECT.
'''
import argparse
import sys

from cyvcf2 import VCF

ANN_FIELDS = ['Allele', 'Annotation', 'Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type', 'Feature_ID',
              'Transcript_BioType', 'Rank', 'HGVS_c', 'HGVS_p', 'cDNA_pos', 'CDS_pos', 'Protein_pos', 'Distance']


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('vcf')
    ap.add_argument('--mane-only', action='store_true')
    ap.add_argument('--impact', default='')
    ap.add_argument('--out', default='parsed_annotations.tsv')
    a = ap.parse_args()
    vcf = VCF(a.vcf)
    tag, fields = None, None
    for h in vcf.header_iter():
        if h['HeaderType'] != 'INFO':
            continue
        if h['ID'] == 'CSQ':
            tag, fields = 'CSQ', h['Description'].split('Format: ')[1].rstrip('"').strip().split('|')
            break
        if h['ID'] == 'ANN' and tag is None:
            tag, fields = 'ANN', ANN_FIELDS
    if tag is None:
        sys.exit('REFUSED: no CSQ (VEP) or ANN (SnpEff) INFO field in the header; annotate first')
    if a.mane_only and (tag != 'CSQ' or 'MANE_SELECT' not in fields):
        sys.exit('REFUSED: --mane-only needs a VEP CSQ with a MANE_SELECT field (run VEP 110+ with --mane_select)')
    impact_key = 'IMPACT' if tag == 'CSQ' else 'Impact'
    keep_impact = {v for v in a.impact.split(',') if v}
    if keep_impact and impact_key not in fields:
        sys.exit(f'REFUSED: the {tag} format has no {impact_key} field')
    n_var = n_rows = 0
    with open(a.out, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\t'.join(['chrom', 'pos', 'ref', 'alt'] + fields) + '\n')
        for v in vcf:
            n_var += 1
            raw = v.INFO.get(tag)
            if not raw:
                continue
            for block in raw.split(','):
                ann = dict(zip(fields, block.split('|')))
                if a.mane_only and not ann.get('MANE_SELECT'):
                    continue
                if keep_impact and ann.get(impact_key) not in keep_impact:
                    continue
                fh.write('\t'.join([v.CHROM, str(v.POS), v.REF, ','.join(v.ALT)] + [ann.get(f, '') for f in fields]) + '\n')
                n_rows += 1
    print(f'{n_var} variants, {n_rows} {tag} blocks written to {a.out}')


if __name__ == '__main__':
    main()
