# Collapse uncertain support
~~~bash
python scripts/edit.py tree.nwk collapsed.nwk collapse --cutoff 70 --zero-floor 0.000001
~~~
Done when the new tree is checked and cutoffs and support measure reported. The source default is 70 for standard bootstrap and 95 for UFBoot2; the last field of dual labels is UFBoot. Zero floor 1e-6 matches the source IQ-TREE example; omit it when inappropriate. Collapsed nodes are soft uncertainty polytomies, not evidence of radiation. Collapsing may change distances; preserve the resolved original. Ladderizing is available via the same script and changes display order only. Arbitrarily resolving polytomies invents topology; see [guidance](../references/caveats.md).
