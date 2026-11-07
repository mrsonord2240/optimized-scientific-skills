# Ligand feature prefilter

```python
import json
from scripts.pharmacophore import feature_family_prefilter
actives = [x.strip() for x in open('actives.smi') if x.strip()]
library = [x.strip() for x in open('library.smi') if x.strip()]
hits = feature_family_prefilter(actives, library, min_shared_fraction=1.0)
json.dump(hits, open('prefilter_hits.json', 'w'), indent=2)
```

The helper is `scripts/pharmacophore.py`. It intersects feature types, not aligned geometry; invalid structures may be skipped. Check rejected known actives and cluster diverse scaffolds. Consensus derivation and an existing RDKit 3D query are described in `references/ligand.md`; feature choices are in `references/selection.md`.

Done when feature prefilter output and a separately validated geometric model are documented.
