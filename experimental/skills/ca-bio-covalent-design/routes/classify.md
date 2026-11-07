# Classify warheads

```python
import json
from scripts.warhead_classifier import classify_warheads
inputs = [line.strip() for line in open('input.smi') if line.strip()]
json.dump([{'smiles': s, 'matches': classify_warheads(s)} for s in inputs], open('warheads.json', 'w'), indent=2)
```

The helper is `scripts/warhead_classifier.py`. Invalid structures return None. SMARTS families overlap; preserve all matches or choose the most specific when a single class is needed. Qualitative tiers are not compound-level rate or selectivity predictions. See `references/warheads.md`.

Done when each input and all recognized matches are recorded.
