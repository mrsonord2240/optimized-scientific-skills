import argparse
from Bio import Phylo

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('--format', required=True, choices=['newick', 'nexus', 'phyloxml', 'nexml'])
args = parser.parse_args()
with open(args.input, encoding='utf-8') as handle:
    for index, tree in enumerate(Phylo.parse(handle, args.format), 1):
        print(f'tree {index}: {len(tree.get_terminals())} tips; rooted={tree.rooted}')
        print('tips:', ', '.join(str(t.name) for t in tree.get_terminals()))
        for node in tree.get_nonterminals():
            print('internal:', 'name=', node.name, 'confidence=', node.confidence, 'comment=', getattr(node, 'comment', None))
