'''SuperSeries/SubSeries check and series-matrix download/parse for one GEO Series.

Inputs:  a GSE accession (relation, matrix) or nothing (selftest).
Outputs: relation -> SuperSeries/SubSeries links from the SOFT family file;
         matrix   -> samples x features shape plus every distinct !Sample_data_processing note;
         selftest -> offline regression test (synthetic non-ASCII fixtures, no network).
Usage:   python geo_series.py relation GSE346738          (streams the SOFT header; nothing saved)
         python geo_series.py matrix GSE470 [--outdir DIR]
         python geo_series.py selftest
Import:  from geo_series import check_super_or_sub_series, parse_series_matrix
Checked: Python 3.12, pandas 3.0, 2026-09-21.

Every gzip.open(..., 'rt') below passes encoding='utf-8', errors='replace': GEO files hold UTF-8
and Python's default on Windows is the locale (cp1252), which raises UnicodeDecodeError.
'''
import argparse
import ast
import gzip
import os
import sys
import tempfile
import urllib.request

import pandas as pd


def parse_series_relation(src):
    '''Read !Series_relation lines from a SOFT family file (.soft.gz): a path or a binary file object.'''
    super_of = []
    sub_of = None
    with gzip.open(src, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('!Series_relation'):
                if 'SuperSeries of' in line:
                    super_of.append(line.split('SuperSeries of: ')[1].strip())
                elif 'SubSeries of' in line:
                    sub_of = line.split('SubSeries of: ')[1].strip()
            if line.startswith(('^PLATFORM', '^SAMPLE')):
                break   # Series block ends here; platform tables can run to millions of lines
    return {'super_of': super_of, 'sub_of': sub_of}


def check_super_or_sub_series(gse):
    # Streamed, not saved: a family file holds every sample (GSE122288's is ~490 MB) and the
    # relation lines sit in the header, so parsing stops at the first ^PLATFORM or ^SAMPLE.
    prefix = gse[:-3] + 'nnn'
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/soft/{gse}_family.soft.gz'
    with urllib.request.urlopen(url, timeout=60) as resp:
        return parse_series_relation(resp)


def download_series_matrix(gse, outdir='.'):
    prefix = gse[:-3] + 'nnn'
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz'
    path = os.path.join(outdir, f'{gse}_matrix.txt.gz')
    urllib.request.urlretrieve(url, path)
    return path


def parse_series_matrix(path):
    metadata = {}
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!'):
                key, *vals = line.rstrip('\n').split('\t')
                metadata[key] = [v.strip('"') for v in vals]
        expr = pd.read_csv(f, sep='\t', index_col=0, comment='!')
    # Series matrix values are whatever submitter chose -- check metadata['!Sample_data_processing']
    return metadata, expr


def selftest():
    '''Offline regression guard for the missing-encoding bug class (no network).'''
    # 1. Static: every gzip.open(..., 'rt') in this file must pass encoding=.
    tree = ast.parse(open(__file__, encoding='utf-8').read())
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == 'open' and getattr(node.func.value, 'id', '') == 'gzip'
                and any(isinstance(a, ast.Constant) and a.value == 'rt' for a in node.args)):
            kw = {k.arg for k in node.keywords}
            assert 'encoding' in kw, f'gzip.open(..., "rt") without encoding= at line {node.lineno}'

    # 2. Dynamic: fixtures with U+0141 (UTF-8 c5 81; 0x81 is undefined in cp1252).
    soft = ('^SERIES = GSE1\n!Series_title = Łódź – café study\n'
            '!Series_relation = SuperSeries of: GSE2\n!Series_relation = SuperSeries of: GSE3\n'
            '^PLATFORM = GPL1\n!Series_relation = SubSeries of: NEVER_READ\n')
    matrix = ('!Series_geo_accession\t"GSE1"\n'
              '!Sample_data_processing\t"Łódź RMA"\t"Łódź RMA"\n'
              '!series_matrix_table_begin\n"ID_REF"\t"GSM1"\t"GSM2"\n"g1"\t1.5\t2.5\n"g2"\t3.5\t4.5\n'
              '!series_matrix_table_end\n')
    with tempfile.TemporaryDirectory() as d:
        sp, mp = os.path.join(d, 't.soft.gz'), os.path.join(d, 't_matrix.txt.gz')
        for p, txt in ((sp, soft), (mp, matrix)):
            with gzip.open(p, 'wb') as f:
                f.write(txt.encode('utf-8'))
            try:   # the fixture must actually trip the bug it guards against
                with gzip.open(p, 'rt', encoding='cp1252') as f:
                    f.read()
                raise AssertionError('fixture does not discriminate: cp1252 read succeeded')
            except UnicodeDecodeError:
                pass
        rel = parse_series_relation(sp)
        assert rel == {'super_of': ['GSE2', 'GSE3'], 'sub_of': None}, rel
        meta, expr = parse_series_matrix(mp)
        assert expr.shape == (2, 2), expr.shape
        assert meta['!Series_geo_accession'] == ['GSE1'], meta
        assert set(meta['!Sample_data_processing']) == {'Łódź RMA'}, meta
    print('selftest OK')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('relation').add_argument('gse')
    p = sub.add_parser('matrix')
    p.add_argument('gse')
    p.add_argument('--outdir', default='.')
    sub.add_parser('selftest')
    a = ap.parse_args()
    if a.cmd == 'selftest':
        selftest()
    elif a.cmd == 'relation':
        print(check_super_or_sub_series(a.gse))
    else:
        meta, expr = parse_series_matrix(download_series_matrix(a.gse, a.outdir))
        print(f'{a.gse}: {expr.shape[0]} features x {expr.shape[1]} samples')
        print('Sample-level data processing notes:')
        for note in sorted(set(meta.get('!Sample_data_processing', []))):
            print(f'  - {note}')


if __name__ == '__main__':
    sys.exit(main())
