# Query similarity search

```python
import json
from scripts.similarity_search import find_similar
query = open('query.smi').readline().strip()
library = [x.strip() for x in open('library.smi') if x.strip()]
json.dump(find_similar(query, library, threshold=0.7, fp_type='ecfp4'), open('hits.json', 'w'), indent=2)
```

The helper is `scripts/similarity_search.py`; it skips invalid library SMILES, so validate and report those inputs first. Threshold 0.7 is a starting heuristic. For Tversky, Dice and alternate vectors see `references/coefficients.md`; confirm substructure-like hits with SMARTS.

Done when ranked hits, fingerprint settings and invalid inputs are recorded.
