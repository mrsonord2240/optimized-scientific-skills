## MAP4 and MHFP6 for Diverse Libraries

For libraries spanning drug-like molecules, natural products, peptides, and metabolites, compare ECFP4 with MAP4 or MHFP6 on task-relevant retrieval benchmarks. MAP4 and MHFP6 use MinHash with atom-pair or circular-substructure shingles, but no universal pairwise-similarity range establishes that ECFP4 is saturated for every mixed library.

### MAP4 installation boundary

Do **not** use `pip install map4`: the historical `map4` name has no supported PyPI distribution in the tested environment, so this Skill does not provide a runnable MAP4 snippet or claim a tested `map4` version. If MAP4 is essential, install and pin a maintained upstream source implementation in a separately reproducible environment; its `tmap` dependency may require a conda/source installation and has no tested Windows build here. Record the source commit, Python/conda platform, and a one-molecule construction check before using its scores. For a separately packaged stereochemistry-aware MinHash option, evaluate `mapchiral` against the task benchmark rather than assuming it is a drop-in MAP4 replacement.

```python
from mhfp.encoder import MHFPEncoder

encoder = MHFPEncoder(2048)
mhfp6 = encoder.encode_mol(mol, radius=3)
```

MHFP6 distance is Jaccard on MinHash, not standard Tanimoto. Use `MHFPEncoder.distance(fp1, fp2)`.
