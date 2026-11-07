import argparse
from itertools import combinations
from Bio import Phylo

p = argparse.ArgumentParser()
p.add_argument('input')
p.add_argument('--taxa', nargs='+', required=True)
a = p.parse_args()
tree = Phylo.read(a.input, 'newick')
tips = tree.get_terminals()
names = [tip.name for tip in tips]
if len(set(names)) != len(names) or not set(a.taxa).issubset(names):
    p.error('unique named tips required; requested taxa must exist')
nodes = [tree.find_any(name=name) for name in a.taxa]
mrca = tree.common_ancestor(*nodes)
print('monophyletic:', bool(tree.is_monophyletic(nodes)))
print('MRCA descendant tips:', sorted(t.name for t in mrca.get_terminals()))
for left, right in combinations(nodes, 2):
    print(left.name, right.name, tree.distance(left, right), sep='\t')
