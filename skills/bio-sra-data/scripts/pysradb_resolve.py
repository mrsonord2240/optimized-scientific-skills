'''Resolve GSE / BioProject / SRX / SRP IDs to SRR run accessions with pysradb.

Inputs:  one or more accessions on the command line (GSE..., PRJNA/PRJEB..., SRP..., SRX...).
Usage:   python scripts/pysradb_resolve.py GSE110009 [PRJEB37378 ...]
         or import: from pysradb_resolve import gse_to_srr, bioproject_to_runs, batch_resolve
Output:  one line per ID with its run count, then run accessions (one per line) on stdout.
'''
# Reference: pysradb 2.2+ (checked 2.5.1) | Verify API if version differs
import sys

import pandas as pd
from pysradb import SRAweb


def gse_to_srr(gse):
    db = SRAweb()
    df = db.gse_to_srp(gse)
    if df.empty:
        return []
    srp = df['study_accession'].iloc[0]
    runs = db.srp_to_srr(srp)
    return runs['run_accession'].tolist()


def bioproject_to_runs(prjna):
    db = SRAweb()
    return db.sra_metadata(prjna, detailed=True)


def batch_resolve(ids):
    db = SRAweb()
    rows = []
    for id in ids:
        try:
            meta = db.sra_metadata(id, detailed=True)
            rows.append(meta)
        except Exception as e:
            print(f'{id}: {e}')
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


if __name__ == '__main__':
    ids = sys.argv[1:] or ['GSE123456']
    for acc in ids:
        if acc.upper().startswith('GSE'):
            srrs = gse_to_srr(acc)
        else:
            srrs = bioproject_to_runs(acc)['run_accession'].tolist()
        print(f'{acc} -> {len(srrs)} SRRs', file=sys.stderr)
        print('\n'.join(srrs))
