# Linker reach prefilter

```python
import json
from scripts.ternary_geometry_screen import screen_linker_library
span = float(open('required_span_angstrom.txt').read().strip())
json.dump(screen_linker_library(span, tolerance_A=2.0), open('linker_reach.json', 'w'), indent=2)
```

The helper is `scripts/ternary_geometry_screen.py`. Measure the required span only after placing the ligand attachment points in a shared structural frame. The standalone-linker conformer range has no protein sterics or interface energy; it is not ternary feasibility. See `references/geometry.md`.

Done when sampled ranges, required span and limitations are recorded.
