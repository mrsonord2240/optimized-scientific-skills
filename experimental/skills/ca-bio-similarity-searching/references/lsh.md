## Large-Library Nearest Neighbor (MHFP6 + LSH Forest)

For large libraries, direct all-pairs comparison becomes expensive. The `mhfp` package provides an LSH-forest helper for approximate nearest-neighbor retrieval over MHFP6 fingerprints. Measure recall and latency against an exact subset for the project dataset.

```python
from mhfp.encoder import MHFPEncoder
from mhfp.lsh_forest import LSHForestHelper

encoder = MHFPEncoder(2048)

def build_index(smiles_list):
    forest = LSHForestHelper()
    fingerprints = []
    for i, smiles in enumerate(smiles_list):
        fp = encoder.encode(smiles, radius=3)
        fingerprints.append(fp)
        forest.add(i, fp)
    forest.index()
    return forest, fingerprints

def query_index(forest, qmol, fingerprints, k=10):
    qfp = encoder.encode_mol(qmol, radius=3)
    return forest.query(qfp, k=k, data=fingerprints)
```

The returned neighbors are approximate in MHFP6 space. Benchmark them against an exact MHFP-distance search on a representative subset before choosing LSH parameters.
