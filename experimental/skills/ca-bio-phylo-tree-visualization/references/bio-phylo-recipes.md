## Bio.Phylo + matplotlib Recipes

Quick text and ASCII inspection, no figure needed:

```python
from Bio import Phylo

tree = Phylo.read('tree.nwk', 'newick')
print(tree)                  # indented text summary
Phylo.draw_ascii(tree)       # ASCII-art diagram, useful in a terminal or log
```

Draw to a vector file (always pass an axes and `do_show=False` for headless/scripted use):

```python
from Bio import Phylo
import matplotlib.pyplot as plt

tree = Phylo.read('tree.nwk', 'newick')
tree.ladderize()             # legibility only; ordering carries NO phylogenetic meaning -- say so in the caption

fig, ax = plt.subplots(figsize=(10, 8))
Phylo.draw(tree, axes=ax, do_show=False)
ax.set_title('Phylogenetic tree (phylogram, branch length = subs/site)')
fig.savefig('tree.pdf', bbox_inches='tight')   # vector: text and lines stay sharp at any size
plt.close(fig)
```

Label tips, and show support with its measure named (never a bare integer):

```python
def tip_only(clade):
    return clade.name if clade.is_terminal() else ''

import re

for clade in tree.get_nonterminals():            # IQ-TREE -B + --alrt: '88.5/92' stays in clade.name, confidence None
    if clade.confidence is None and clade.name and re.fullmatch(r'[\d.]+/[\d.]+', clade.name):
        clade.sh_alrt, clade.ufboot = (float(v) for v in clade.name.split('/'))   # SH-aLRT / UFBoot (last)
        clade.name = None                        # otherwise the default label_func prints the raw string

def support_label(clade):
    # the measure MUST be stated in the legend/caption
    if clade.is_terminal():
        return ''
    if getattr(clade, 'ufboot', None) is not None:
        return f'{clade.sh_alrt:.0f}/{clade.ufboot:.0f}'
    return f'{clade.confidence:.0f}' if clade.confidence is not None else ''

if not any(support_label(c) for c in tree.get_nonterminals()):
    print('WARNING: no internal clade has readable support -- the figure will show none')
fig, ax = plt.subplots(figsize=(12, 10))
Phylo.draw(tree, axes=ax, do_show=False, label_func=tip_only, branch_labels=support_label)
ax.set_title('Node support: SH-aLRT (%) / UFBoot (%)')   # name the measure(s) the file actually holds
fig.savefig('supported_tree.svg', bbox_inches='tight')
plt.close(fig)
```

Color branches by group (set `.color` on the MRCA clade; descendants inherit it). ROOT FIRST: `common_ancestor` on an arbitrarily-rooted Newick read (the common case -- IQ-TREE etc. write an unrooted trifurcation) can return a basal node spanning nearly the whole tree instead of the intended clade, so root on an outgroup (tree-manipulation) before calling it, and check the returned tip set:

```python
import re

tree.root_with_outgroup({'name': 'OutA'}, {'name': 'OutB'})   # root BEFORE any common_ancestor call -- see tree-manipulation

intended_tips = {'Homo_sapiens', 'Pan_troglodytes', 'Gorilla_gorilla', 'Pongo_abelii'}   # the clade you actually mean
mrca = tree.common_ancestor({'name': 'Homo_sapiens'}, {'name': 'Pongo_abelii'})
mrca_tips = set(t.name for t in mrca.get_terminals())
if mrca_tips != intended_tips:                 # FAIL LOUDLY -- do not silently color the wrong clade
    raise ValueError(f'MRCA gave {sorted(mrca_tips)}, not the intended clade {sorted(intended_tips)} '
                      f'-- check the outgroup/rooting and the tip names passed to common_ancestor')
mrca.color = 'red'

for clade in tree.get_nonterminals():          # clear raw support strings so they are not drawn as node labels
    if clade.name and re.fullmatch(r'[\d.]+/[\d.]+', clade.name):
        clade.name = None

fig, ax = plt.subplots(figsize=(10, 8))
Phylo.draw(tree, axes=ax, do_show=False)       # as_phyloxml() is needed only to EXPORT colors to phyloXML
fig.savefig('colored_tree.pdf', bbox_inches='tight')
plt.close(fig)
```

Scale the panel to tip count so labels stay legible, keeping the branch-length axis (Bio.Phylo has no scale-bar artist, so the x axis IS the scale):

```python
n_tips = len(tree.get_terminals())
MAX_TIPS = 200                                  # Bio.Phylo/matplotlib ceiling: no radial/ring layout, so labels stay illegible past this at any panel height (320 tips checked)
if n_tips > MAX_TIPS:
    raise ValueError(f'{n_tips} tips > {MAX_TIPS}: Bio.Phylo cannot draw legible labels -- use ggtree circular/fan '
                     f'or iTOL (raise MAX_TIPS only to accept an illegible figure)')
if n_tips > 150:
    print('>~150 tips: consider a circular layout or strips/rings (ggtree, iTOL) instead of a taller panel')
height = min(max(8, n_tips * 0.25), 40)         # ~0.25 in/tip keeps ~6-8 pt labels from colliding; capped

fig, ax = plt.subplots(figsize=(10, height))
Phylo.draw(tree, axes=ax, do_show=False)
ax.set_yticks([])                               # hide only the meaningless y ticks and spines
ax.spines[['left', 'top', 'right']].set_visible(False)
fig.savefig('scaled_tree.pdf', bbox_inches='tight')
plt.close(fig)
```

For circular/fan/unrooted layouts, metadata heatmaps, dual support, or BEAST HPD bars, Bio.Phylo cannot help -- route to ggtree + treeio (R), ETE4, or iTOL.
