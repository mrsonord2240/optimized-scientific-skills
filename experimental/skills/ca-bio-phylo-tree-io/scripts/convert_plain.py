import argparse
from pathlib import Path
from Bio import Phylo

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('output')
parser.add_argument('--from-format', required=True, choices=['newick', 'nexus', 'phyloxml', 'nexml'])
parser.add_argument('--to-format', required=True, choices=['newick', 'nexus', 'phyloxml', 'nexml'])
args = parser.parse_args()
if Path(args.input).resolve() == Path(args.output).resolve():
    parser.error('input and output must be different files')
with open(args.input, encoding='utf-8') as source, open(args.output, 'x', encoding='utf-8') as dest:
    n = Phylo.write(Phylo.parse(source, args.from_format), dest, args.to_format)
print(f'wrote {n} tree(s); verify labels and metadata')
