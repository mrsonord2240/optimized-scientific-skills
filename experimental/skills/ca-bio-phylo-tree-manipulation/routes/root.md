# Root an existing tree
~~~bash
python scripts/edit.py tree.nwk rooted.nwk outgroup --ingroup I1 I2 I3 --outgroup OutA OutB
~~~
Done when the new tree is written, separately declared ingroup is checked and root uncertainty is recorded. Outgroups should be close and monophyletic; the script roots temporarily on an ingroup tip and splits the outgroup stem. Check support against bipartitions after rerooting. For a clock-limited fallback use the midpoint operation; long branches can mislead. Compare methods, including MAD, MinVar and likelihood alternatives in [rooting guidance](../references/rooting.md). Agreement alone is not root confidence.
