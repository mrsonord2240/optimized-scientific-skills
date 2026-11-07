# Subset and inspect a tree
~~~bash
python scripts/inspect.py tree.nwk --taxa Human Mouse Rat
python scripts/edit.py tree.nwk subset.nwk prune --keep Human Mouse Rat
python scripts/inspect.py subset.nwk --taxa Human Mouse Rat
~~~
Done when requested tips remain and known surviving distances agree before and after pruning. Input names must be unique. MRCA descendants can include extra tips; non-monophyletic targets yield an induced subtree, not a clade. Retained support is not subset support. See [caveats](../references/caveats.md).
