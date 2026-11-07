#!/usr/bin/env python3
'''
SpliceVault (300K-RNA) lookup of the most likely mis-splicing events for a splice-site SNV, from the
public Ensembl-hosted table (the same file the SpliceVault VEP plugin reads). GRCh38 only here.

Usage: python splicevault_lookup.py chr17 7674292 T C [--transcript ENST00000269305]
Research-use decision support: empirical events seen in reference RNA-seq, not a diagnosis.
Only the variant coordinates are sent (HTTP range requests to ftp.ensembl.org); nothing else leaves the machine.
Checked with pysam 0.24.1 on SpliceVault_data_GRCh38.tsv.gz (Ensembl current_variation, read 2026-09-20).
'''
import argparse
import pandas as pd
import pysam

BASE = 'https://ftp.ensembl.org/pub/current_variation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz'
COLS = ['chrom', 'pos', 'ref', 'alt', 'transcript', 'site_type', 'site', 'spliceai_delta',
        'out_of_frame', 'top_events', 'site_sample_count', 'site_max_depth']


def parse_events(text):
    '''"Top1:CA;+47;0.4%;Frameshift|Top2:ES;3;0.03%;inFrame" -> rank, type (ES exon skip, CA/CD cryptic
    acceptor/donor), impact (skipped exon(s) or nt from the annotated site), % supporting samples, frame.'''
    rows = []
    for ev in text.split('|'):
        rank, rest = ev.split(':', 1)
        typ, impact, pct, frame = rest.split(';')
        rows.append({'rank': rank, 'type': typ, 'impact': impact, 'pct_samples': pct, 'frame': frame})
    return pd.DataFrame(rows)


def splicevault(chrom, pos, ref, alt, path=BASE, index=None):
    '''All SpliceVault rows (one per transcript) for chrom:pos ref>alt. Empty frame = variant not at a
    catalogued splice site (SpliceVault only covers SNVs at/near annotated donors and acceptors).
    `path` may be a local copy of the .tsv.gz; `index` the local .tbi (needed when path is a URL).'''
    tbx = pysam.TabixFile(path, index=index) if index else pysam.TabixFile(path)
    chrom = chrom if chrom.startswith('chr') else 'chr' + chrom
    rows = [r.split('\t') for r in tbx.fetch(chrom, int(pos) - 1, int(pos))]
    df = pd.DataFrame(rows, columns=COLS)
    return df[(df['ref'] == ref) & (df['alt'] == alt)].reset_index(drop=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('chrom'); ap.add_argument('pos', type=int); ap.add_argument('ref'); ap.add_argument('alt')
    ap.add_argument('--transcript', help='keep one transcript (use the MANE Select / canonical ID)')
    ap.add_argument('--tsv', default=BASE); ap.add_argument('--tbi', help='local .tbi when --tsv is a URL')
    a = ap.parse_args()
    df = splicevault(a.chrom, a.pos, a.ref, a.alt, a.tsv, a.tbi)
    if a.transcript:
        df = df[df['transcript'].str.startswith(a.transcript)]
    if df.empty:
        print('no SpliceVault entry (not a catalogued splice-site SNV, or wrong REF/ALT/build)')
    for r in df.itertuples():
        print(f'{r.transcript} {r.site_type} at {r.site}; SpliceAI delta {r.spliceai_delta}; '
              f'out-of-frame {r.out_of_frame}; {r.site_sample_count} samples')
        print(parse_events(r.top_events).to_string(index=False))
