'''Parsimony grouping + PICKED protein-group FDR on a labeled peptide->protein map.

Self-contained demo (no input files): a small INLINE peptide-to-protein mapping
including decoy proteins, so the picked-FDR logic runs end to end. With a real
search, replace SAMPLE_MAP with the peptide->protein lists parsed from an idXML
(pyopenms.IdXMLFile) or a MaxQuant proteinGroups.txt, replace the toy scores
with the inference engine's group probabilities, and set DECOY_PREFIX to the
tool's prefix (MaxQuant: 'REV__').

Indistinguishable proteins (identical observed peptides) are collapsed into one
node BEFORE the greedy loop, so they are reported as one group whatever order the
search listed them in, and are never mislabelled as subsumable.

Expected output: 7 of 8 candidate groups kept, P_C subsumable; P_E and P_E-2
reported as one group led by P_E; DECOY_P_D wins its pick against P_D, so P_B
(ranked after it) gets q = 0.33 and 2 target groups pass. Three decoy groups are
far too few for a real FDR estimate, so the demo prints a warning.
'''
# Reference: pandas 2.2+ | Verify API if version differs
import pandas as pd

DECOY_PREFIX = 'DECOY_'  # marks reversed/shuffled decoy proteins
FDR_CUTOFF = 0.01  # 1% protein-group FDR, the community standard for the LIST
MIN_DECOYS = 10  # below this the decoy count cannot support an FDR estimate

# peptide -> proteins it maps to. 'unique' peptides map to one protein, 'shared'
# (degenerate) peptides map to several. Decoy proteins carry the DECOY_ prefix.
SAMPLE_MAP = {
    'AAAEUNIQUER': ['P_A'],
    'BBBAUNIQUEK': ['P_A'],
    'SHAREDABK': ['P_A', 'P_B'],
    'CCCBUNIQUEK': ['P_B'],
    'DDDSUBSETK': ['P_C', 'P_A'],          # P_C's only peptide is also P_A's -> subsumable
    'EEEONLYK': ['P_D'],
    'FFFISOFORMK': ['P_E-2', 'P_E'],       # P_E and its isoform P_E-2 share every observed peptide
    'GGGISOFORMR': ['P_E-2', 'P_E'],       # -> indistinguishable, isoform listed first on purpose
    'HHHISOFORMK': ['P_E-2', 'P_E'],
    'ZZZDECOYK': [DECOY_PREFIX + 'P_A'],   # decoy counterpart of target P_A; picking pairs them
    'YYYDECOYR': [DECOY_PREFIX + 'P_A'],
    'XXXDECOYK': [DECOY_PREFIX + 'P_B'],   # decoy counterpart of target P_B
    'WWWDECOYK': [DECOY_PREFIX + 'P_D'],   # two decoy peptides beat P_D's one: the decoy wins the pick
    'VVVDECOYR': [DECOY_PREFIX + 'P_D'],
}


def protein_nodes(peptide_protein_map):
    '''Collapse proteins with identical observed-peptide sets into one node (indistinguishable group).'''
    protein_peptides = {}
    for pep, prots in peptide_protein_map.items():
        for p in prots:
            protein_peptides.setdefault(p, set()).add(pep)
    nodes = {}
    for prot, peps in protein_peptides.items():
        nodes.setdefault(frozenset(peps), []).append(prot)
    # canonical accession (no -N isoform suffix) first, then alphabetical: independent of input order
    return {peps: sorted(prots, key=lambda a: ('-' in a, a)) for peps, prots in nodes.items()}


def apply_parsimony(peptide_protein_map):
    nodes = protein_nodes(peptide_protein_map)
    all_peptides = set(peptide_protein_map)
    covered = set()
    selected = []
    while covered != all_peptides:
        # most new peptides; ties broken by total peptides, then accession (deterministic)
        best = max(nodes, key=lambda n: (len(n - covered), len(n), [-ord(c) for c in nodes[n][0]]))
        if not best - covered:
            break
        selected.append(best)
        covered |= best
    dropped = [p for n, prots in nodes.items() if n not in selected for p in prots]
    return selected, nodes, sorted(dropped)


def build_groups(selected, nodes, peptide_protein_map):
    groups = []
    for peptides in selected:
        prots = nodes[peptides]
        n_unique = sum(1 for pep in peptides if set(peptide_protein_map[pep]) <= set(prots))
        groups.append({
            'leading_protein': prots[0],
            'accessions': prots,
            'n_peptides': len(peptides),
            'n_unique_peptides': n_unique,   # unique to the GROUP (MaxQuant's meaning)
            'score': len(peptides) + n_unique,  # toy score; use engine probability with real data
            'is_decoy': all(a.startswith(DECOY_PREFIX) for a in prots),
        })
    return groups


def picked_group_fdr(groups, decoy_prefix, min_decoys=MIN_DECOYS):
    # PICKED FDR: pair each target group with its decoy counterpart (same accessions
    # minus the prefix) and keep only the higher-scoring of the pair before counting.
    n_decoy = sum(g['is_decoy'] for g in groups)
    if n_decoy == 0:
        raise ValueError(f'no decoy groups with prefix {decoy_prefix!r}; keep decoys upstream or fix the prefix')
    if n_decoy < min_decoys:
        print(f'WARNING: only {n_decoy} decoy groups; the protein FDR estimate is not meaningful')
    by_base = {}
    for g in groups:
        base = frozenset(a.replace(decoy_prefix, '') for a in g['accessions'])
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
    for g in reversed(picked):  # enforce monotone q-values from the bottom up
        running_min = min(running_min, g['fdr'])
        g['qvalue'] = running_min
    return picked


selected, nodes, dropped = apply_parsimony(SAMPLE_MAP)
print(f'Parsimony kept {len(selected)} groups out of {len(nodes)} candidate groups')
print(f'Dropped (subsumable): {dropped}')

groups = build_groups(selected, nodes, SAMPLE_MAP)
picked = picked_group_fdr(groups, DECOY_PREFIX)

report = pd.DataFrame(picked).sort_values('score', ascending=False)
passing = report[(~report['is_decoy']) & (report['qvalue'] <= FDR_CUTOFF)]
print(f'\nProtein groups passing {FDR_CUTOFF:.0%} picked-group FDR: {len(passing)}')
print(report[['leading_protein', 'accessions', 'n_peptides', 'n_unique_peptides', 'is_decoy', 'qvalue']].to_string(index=False))
