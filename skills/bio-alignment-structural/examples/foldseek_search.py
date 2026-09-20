'''Search a query structure against AlphaFoldDB or a custom Foldseek database.

Returns best hits with TM-score, e-value, and structural alignment metadata.
--alignment-type 2 (3Di+AA Gotoh local) is Foldseek's default and is the right
choice for the same-fold regime; --alignment-type 1 (TMalign) refines top hits
with a full global TM-score at higher cost. E-values are not meaningful under type 1
(true homologs report 0.88-0.99), so confident_hits() filters on alntmscore there.
'''
# Reference: foldseek 8+ (checked on 10.941cd33) | Verify CLI flags if version differs

import subprocess
import csv
import math

def foldseek_search(query_pdb, database, output_m8, tmp_dir='tmp/', alignment_type=2, max_seqs=200):
    cmd = [
        'foldseek', 'easy-search', query_pdb, database, output_m8, tmp_dir,
        '--alignment-type', str(alignment_type),
        '--max-seqs', str(max_seqs),
        '--format-output', 'query,target,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,alnlen,pident',
    ]
    subprocess.run(cmd, check=True)

def parse_results(output_m8):
    columns = ['query', 'target', 'evalue', 'bits', 'alntmscore', 'qtmscore', 'ttmscore', 'lddt', 'alnlen', 'pident']
    hits = []
    with open(output_m8) as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) != len(columns):
                print(f'Warning: skipping row with {len(row)} fields (expected {len(columns)})')
                continue
            hit = dict(zip(columns, row))
            for col in ['evalue', 'bits', 'alntmscore', 'qtmscore', 'ttmscore', 'lddt', 'pident']:
                value = hit[col]
                hit[col] = float(value) if value not in ('', 'NA', 'nan') else math.nan
            hit['alnlen'] = int(hit['alnlen']) if hit['alnlen'].isdigit() else 0
            hits.append(hit)
    return hits

def confident_hits(hits, alignment_type=2):
    '''alnTM > 0.5 always; also E-value < 1e-3 for alignment type 2 (E-value is uninformative under type 1).'''
    keep = [h for h in hits if not math.isnan(h['alntmscore']) and h['alntmscore'] > 0.5]
    if alignment_type != 1:
        keep = [h for h in keep if h['evalue'] < 1e-3]
    return keep

if __name__ == '__main__':
    alignment_type, max_seqs = 2, 200
    foldseek_search('query.pdb', '/path/to/afdb', 'result.m8', alignment_type=alignment_type, max_seqs=max_seqs)
    hits = parse_results('result.m8')
    if len(hits) >= max_seqs:
        print(f'Warning: {len(hits)} rows = the --max-seqs cap; the database has more hits. Raise max_seqs to count them all.')

    print('Top 20 Foldseek hits:')
    print(f'{"target":<25} {"evalue":>10} {"alnTM":>7} {"qTM":>7} {"lddt":>6} {"%id":>5}')
    for hit in hits[:20]:
        print(f'{hit["target"][:25]:<25} {hit["evalue"]:>10.2e} {hit["alntmscore"]:>7.3f} '
              f'{hit["qtmscore"]:>7.3f} {hit["lddt"]:>6.3f} {hit["pident"]:>5.1f}')

    confident_homologs = confident_hits(hits, alignment_type)
    criterion = 'alnTM > 0.5' if alignment_type == 1 else 'alnTM > 0.5, e-value < 1e-3'
    print(f'\nConfident structural hits ({criterion}) among {len(hits)} rows: {len(confident_homologs)}')
