"""q-values from any engine's PSM table, for a single CONCATENATED target-decoy search.

Input : tab-separated table with one column each for scan, score and protein.
        score must be HIGHER-is-better: Sage sage_discriminant_score or Comet xcorr as is;
        E-values (Comet e-value, MS-GF+ SpecEValue) as -log10(E-value).
        Not for separate target/decoy searches: use examples/separate_search_fdr.py.
Output: target PSMs at q <= 0.01 (--out), and the counts on stdout.
Usage : python table_fdr.py search_results.tsv --scan scan --score score --protein protein [--out kept.tsv]
Checked on pandas 3.0.5, numpy 2.5.3.
"""
import argparse

import pandas as pd

DECOY_PREFIXES = ('decoy_', 'rev_', 'xxx_')   # compared lower-cased: DECOY_, rev_ (Sage), XXX_

ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
ap.add_argument('table')
ap.add_argument('--scan', default='scan', help='column naming the spectrum')
ap.add_argument('--score', default='score', help='higher-is-better score column')
ap.add_argument('--protein', default='protein', help='column with the protein accession(s)')
ap.add_argument('--out', help='write the 1%% list here (TSV)')
args = ap.parse_args()

psms = pd.read_csv(args.table, sep='\t').rename(
    columns={args.scan: 'scan', args.score: 'score', args.protein: 'protein'})
psms['is_decoy'] = psms['protein'].str.lower().str.startswith(DECOY_PREFIXES)
if not psms['is_decoy'].any():
    raise ValueError(f'no decoy PSMs recognised in {len(psms)} rows: check the decoy prefix, or the table was already decoy-filtered '
                     f'(with no decoys the smallest reachable q is 1/{len(psms)} = {1/len(psms):.3f})')
# one best hit per spectrum (Comet .txt writes 5 rows per scan by default)
psms = psms.sort_values('score', ascending=False).drop_duplicates('scan').reset_index(drop=True)

# concatenated target-decoy competition: each decoy above threshold estimates one false target
targets = (~psms['is_decoy']).cumsum()
decoys = psms['is_decoy'].cumsum()
psms['fdr'] = (decoys + 1) / targets.clip(lower=1)   # +1: zero decoys is not zero FDR (OpenMS conservative default)
psms['qvalue'] = psms['fdr'][::-1].cummin()[::-1]   # running min from the bottom -> monotone q-values

kept = psms[(psms['qvalue'] <= 0.01) & (~psms['is_decoy'])]   # 1% list-level FDR
print(f'{len(psms)} spectra, {int(psms["is_decoy"].sum())} decoy hits, {len(kept)} target PSMs at q <= 0.01')
if args.out:
    kept.to_csv(args.out, sep='\t', index=False)
