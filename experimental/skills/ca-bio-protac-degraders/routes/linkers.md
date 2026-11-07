# Enumerate connected candidates

```python
import json
from scripts.protac_enumerate import enumerate_linkers
target = open('target_fragment.smi').readline().strip()
e3 = open('e3_fragment.smi').readline().strip()
json.dump(enumerate_linkers(target, e3), open('candidates.json', 'w'), indent=2)
```

The helper is `scripts/protac_enumerate.py`. Each fragment needs exactly one single-bond attachment dummy; linkers require mapped endpoints. This builds connected structures, not validated synthetic routes or ternary predictions. For E3 and linker decisions see `references/selection.md` and `references/geometry.md`.

Done when candidate connectivity and assumptions are recorded for review.
