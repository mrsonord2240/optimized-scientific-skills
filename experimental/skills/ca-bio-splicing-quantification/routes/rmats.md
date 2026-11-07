# rMATS event PSI

```python
from scripts.quantify_splicing import parse_rmats_output
parse_rmats_output('rmats_output', 'SE', min_junction_reads=20).to_csv('reliable_se.csv', index=False)
```

Use an existing rMATS output directory; for BAM quantification and correct library flags see `references/rmats.md`. SE, A5SS, A3SS, MXE and RI are supported, with JC or JCEC counts. Each replicate in both groups must meet the selected minimum. Use per-event effective lengths and preserve group sign; see `references/selection.md` and `references/interpretation.md`.

Done when per-replicate PSI and requested filtered events are reported.
