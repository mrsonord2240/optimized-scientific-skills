import argparse
import json
import dendropy

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('output')
parser.add_argument('--format', required=True, choices=['newick', 'nexus'])
args = parser.parse_args()
with open(args.input, encoding='utf-8') as handle:
    trees = dendropy.TreeList.get(file=handle, schema=args.format, extract_comment_metadata=True)
records = []
for index, tree in enumerate(trees, 1):
    for node in tree.preorder_node_iter():
        tips = sorted(leaf.taxon.label for leaf in node.leaf_iter())
        annotations = {annotation.name: annotation.value for annotation in node.annotations}
        records.append({'tree': index, 'tips': tips, 'annotations': annotations})
with open(args.output, 'w', encoding='utf-8') as handle:
    json.dump(records, handle, indent=2, ensure_ascii=False, default=str)
    handle.write('\n')
