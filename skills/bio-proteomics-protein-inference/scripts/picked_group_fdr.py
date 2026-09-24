'''Picked protein-group FDR on a table of groups (target/decoy pairing by exact accession set).

Purpose: pair each target group with its decoy counterpart (same accessions minus the decoy prefix),
keep the higher-scoring member of each pair, rank the picked set, count decoys, and return the target
groups at q <= 0.01. For idXML prefer the built-in FalseDiscoveryRate().applyPickedProteinFDR (see
references/pyopenms-basic-inference.md); for large data see the kusterlab `picked_group_fdr` package.
Import: `from picked_group_fdr import picked_group_fdr` (groups = dicts with 'accessions' (list of str),
'score' (higher = better), 'is_decoy').
Inputs (CLI): a TSV with an `accessions` column (;-joined) and a score column (higher = better), e.g.
the --out of epifany_inference.py; the decoy prefix is required (DECOY_ OpenMS/Comet, rev_ Philosopher, REV__ MaxQuant).
Usage: python picked_group_fdr.py groups.tsv --decoy-prefix DECOY_ [--score-col probability] [--out passing.tsv]
'''


def picked_group_fdr(groups, decoy_prefix, min_decoys=10):
    # groups: list of dicts with 'accessions' (str), 'score' (higher = better), 'is_decoy'
    n_decoy = sum(g['is_decoy'] for g in groups)
    if n_decoy == 0 or not any(a.startswith(decoy_prefix) for g in groups for a in g['accessions']):
        raise ValueError(f'no decoy groups with prefix {decoy_prefix!r}; keep decoys upstream or fix the prefix')
    if n_decoy < min_decoys:
        print(f'WARNING: only {n_decoy} decoy groups; the protein FDR estimate is not meaningful')
    by_base = {}
    for g in groups:
        base = frozenset(a.replace(decoy_prefix, '') for a in g['accessions'])
        # keep only the higher-scoring of the target/decoy pair (the 'pick')
        if base not in by_base or g['score'] > by_base[base]['score']:
            by_base[base] = g
    picked = sorted(by_base.values(), key=lambda g: g['score'], reverse=True)

    targets = decoys = 0
    for g in picked:
        if g['is_decoy']:
            decoys += 1
        else:
            targets += 1
        g['fdr'] = decoys / targets if targets else 1.0
    running_min = 1.0
    for g in reversed(picked):  # monotone q-values from the bottom up
        running_min = min(running_min, g['fdr'])
        g['qvalue'] = running_min
    return [g for g in picked if not g['is_decoy'] and g['qvalue'] <= 0.01]


if __name__ == '__main__':
    import argparse
    import csv

    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('tsv')
    ap.add_argument('--decoy-prefix', required=True)
    ap.add_argument('--score-col', default='probability')
    ap.add_argument('--out')
    args = ap.parse_args()
    with open(args.tsv, encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    groups = []
    for r in rows:
        accs = r['accessions'].split(';')
        groups.append({'accessions': accs, 'score': float(r[args.score_col]),
                       'is_decoy': all(a.startswith(args.decoy_prefix) for a in accs)})
    passing = picked_group_fdr(groups, args.decoy_prefix)
    print(len(passing), 'target groups at 1% picked-group FDR')
    if args.out:
        with open(args.out, 'w', encoding='utf-8') as fh:
            fh.write('accessions\tscore\tqvalue\n')
            for g in passing:
                fh.write(';'.join(g['accessions']) + f"\t{g['score']}\t{g['qvalue']}\n")
