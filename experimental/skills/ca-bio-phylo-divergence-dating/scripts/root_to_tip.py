'''Root-to-tip regression: a scripted TempEst replacement for the temporal-signal check.

TempEst (Rambaut et al. 2016) is a GUI with no CLI, so this regresses each tip's
root-to-tip genetic distance against its sampling date directly from an IQ-TREE
Newick file, trying every edge as the root and keeping the rooting that maximizes
R^2 -- the same heuristic TempEst's "best-fitting root" search uses. A genuine
temporal signal gives a positive slope (the slope estimates the substitution
rate); R^2 is exploratory only (tips are not independent data points) -- treat
near-zero / much below ~0.2 as a red flag, not a formal pass. This regression is
necessary but not sufficient: always follow it with date_randomization.py before
trusting a tip-dated result.

Usage:
  python root_to_tip.py TREEFILE DATES.TSV
    TREEFILE   Newick tree (e.g. from `iqtree2 -s seqs.fa -m GTR+G --prefix rttree`)
    DATES.TSV  tip <TAB> decimal-year-date, one row per tip (same format --date takes)
'''
# Reference: DendroPy 5.0, SciPy 1.18 | Verify import paths if version differs

import sys

import dendropy
from scipy import stats


def load_dates(path):
    dates = {}
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            tip, date = line.split('\t')
            dates[tip] = float(date)
    return dates


def root_to_tip_r2(tree, dates):
    '''Regress each tip's root-to-tip patristic distance against its sampling date.'''
    tree.calc_node_root_distances(return_leaf_distances_only=True)
    xs, ys = [], []
    for leaf in tree.leaf_node_iter():
        name = leaf.taxon.label.replace(' ', '_')
        if name not in dates:
            continue
        xs.append(dates[name])
        ys.append(leaf.root_distance)
    if len(xs) < 3:
        raise ValueError(f'only {len(xs)} tips matched dates file; need >= 3')
    slope, intercept, r, _, _ = stats.linregress(xs, ys)
    return slope, intercept, r ** 2, len(xs)


def best_rooting(treefile, dates):
    '''Try every edge as the root; keep the rooting with the highest R^2
    (TempEst's heuristic-residual-mean-squared root search).'''
    tns = dendropy.TaxonNamespace()
    probe = dendropy.Tree.get(path=treefile, schema='newick', taxon_namespace=tns)
    n_edges = sum(1 for e in probe.preorder_edge_iter() if e.head_node.parent_node is not None)
    best = None
    for i in range(n_edges):
        tree = dendropy.Tree.get(path=treefile, schema='newick', taxon_namespace=tns)
        edges = [e for e in tree.preorder_edge_iter() if e.head_node.parent_node is not None]
        edge = edges[i]
        half = (edge.length or 0.0) / 2.0
        tree.reroot_at_edge(edge, length1=half, length2=half, update_bipartitions=False)
        result = root_to_tip_r2(tree, dates)
        if best is None or result[2] > best[2]:
            best = result
    return best


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    treefile, datefile = sys.argv[1], sys.argv[2]
    dates = load_dates(datefile)
    slope, intercept, r2, n = best_rooting(treefile, dates)
    x_intercept = -intercept / slope if slope else float('nan')
    print(f'n_tips_used={n}')
    print(f'slope (rate, subs/site/yr) = {slope:.6e}')
    print(f'intercept = {intercept:.6f}')
    print(f'R^2 = {r2:.4f}')
    print(f'x-intercept (TMRCA estimate) = {x_intercept:.3f}')
    if slope <= 0:
        print('WARNING: slope is not positive -- no usable temporal signal or wrong root.')
    elif r2 < 0.2:
        print('WARNING: R^2 << 0.2 -- weak signal; still run date_randomization.py before trusting a date.')
