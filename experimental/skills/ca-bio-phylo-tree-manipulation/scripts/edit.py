import argparse
from pathlib import Path
from Bio import Phylo

p = argparse.ArgumentParser()
p.add_argument('input')
p.add_argument('output')
p.add_argument('operation', choices=['midpoint', 'outgroup', 'prune', 'collapse', 'ladderize'])
p.add_argument('--keep', nargs='+')
p.add_argument('--ingroup', nargs='+')
p.add_argument('--outgroup', nargs='+')
p.add_argument('--cutoff', type=float, default=70)
p.add_argument('--zero-floor', type=float)
a = p.parse_args()
if Path(a.input).resolve() == Path(a.output).resolve():
    p.error('input and output must differ')
tree = Phylo.read(a.input, 'newick')
if a.operation == 'midpoint':
    tree.root_at_midpoint()
elif a.operation == 'outgroup':
    if not a.ingroup or not a.outgroup or set(a.ingroup) & set(a.outgroup):
        p.error('declare disjoint --ingroup and --outgroup names')
    tips = tree.get_terminals()
    names = [t.name for t in tips]
    if len(set(names)) != len(names) or not set(a.ingroup + a.outgroup).issubset(names):
        p.error('unique named tips required; declared taxa must exist')
    tree.root_with_outgroup({'name': a.ingroup[0]})
    outgroups = [tree.find_any(name=n) for n in a.outgroup]
    if not tree.is_monophyletic(outgroups):
        p.error('outgroup is not monophyletic after temporary ingroup root')
    stem = tree.common_ancestor(*outgroups).branch_length
    if stem is None:
        p.error('outgroup stem has no branch length')
    tree.root_with_outgroup(*outgroups, outgroup_branch_length=stem / 2)
    if not tree.is_monophyletic([tree.find_any(name=n) for n in a.ingroup]):
        print('WARNING: declared ingroup is not monophyletic; reconsider rooting')
elif a.operation == 'prune':
    if not a.keep or len(a.keep) < 2:
        p.error('prune requires at least two --keep names')
    names = [t.name for t in tree.get_terminals()]
    if len(set(names)) != len(names) or not set(a.keep).issubset(names):
        p.error('tip names must be unique and requested tips must exist')
    for tip in list(tree.get_terminals()):
        if tip.name not in a.keep:
            tree.prune(tip)
elif a.operation == 'collapse':
    def support(c):
        if c.confidence is not None:
            return c.confidence
        if c.name and '/' in c.name:
            return float(c.name.split('/')[-1])
        return None
    if not any(support(c) is not None for c in tree.get_nonterminals()):
        p.error('no readable support')
    tree.collapse_all(lambda c: support(c) is not None and support(c) < a.cutoff)
    if a.zero_floor is not None:
        tree.collapse_all(lambda c: c.branch_length is not None and c.branch_length <= a.zero_floor)
else:
    tree.ladderize()
with open(a.output, 'x', encoding='utf-8') as dest:
    Phylo.write(tree, dest, 'newick')
