import argparse
import re
from pathlib import Path
from Bio import Phylo

p = argparse.ArgumentParser()
p.add_argument('input')
p.add_argument('--output')
p.add_argument('--ascii', action='store_true')
p.add_argument('--title', default='Phylogram')
p.add_argument('--support-label', help='Name the support measure and order')
p.add_argument('--ladderize', action='store_true')
p.add_argument('--clade', nargs='+', help='Complete intended descendant tip set on a previously rooted tree')
a = p.parse_args()
if not a.ascii and not a.output:
    p.error('provide --output or --ascii')
if a.output and Path(a.input).resolve() == Path(a.output).resolve():
    p.error('input and output must differ')
tree = Phylo.read(a.input, 'newick')
if a.ladderize:
    tree.ladderize()
if a.clade:
    tips = tree.get_terminals()
    names = [tip.name for tip in tips]
    if len(names) != len(set(names)) or not set(a.clade).issubset(names):
        p.error('unique named tips required; requested clade tips must exist')
    mrca = tree.common_ancestor(*[tree.find_any(name=n) for n in a.clade])
    if {t.name for t in mrca.get_terminals()} != set(a.clade):
        p.error('MRCA has additional tips; check rooting and intended clade')
    mrca.color = 'red'
if a.ascii:
    Phylo.draw_ascii(tree)
if a.output:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = len(tree.get_terminals())
    if n > 200:
        p.error('more than 200 tips: use circular/fan or alternative tool')
    for c in tree.get_nonterminals():
        if c.confidence is None and c.name and re.fullmatch(r'[0-9.]+/[0-9.]+', c.name):
            c.sh_alrt, c.ufboot = (float(v) for v in c.name.split('/'))
            c.name = None
    def support(c):
        if c.is_terminal():
            return ''
        if getattr(c, 'ufboot', None) is not None:
            return f'{c.sh_alrt:.0f}/{c.ufboot:.0f}'
        return f'{c.confidence:.0f}' if c.confidence is not None else ''
    fig, ax = plt.subplots(figsize=(10, min(max(8, n * .25), 40)))
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: c.name if c.is_terminal() else '', branch_labels=support if a.support_label else None)
    ax.set_title(a.title + ((' - ' + a.support_label) if a.support_label else ''))
    ax.set_xlabel('Branch length (state units in caption)')
    with open(a.output, 'xb') as dest:
        fig.savefig(dest, format=Path(a.output).suffix.lstrip('.'), bbox_inches='tight')
    plt.close(fig)
