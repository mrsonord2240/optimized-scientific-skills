# Distance correction models

Read when converting alignment identity into evolutionary distances for publication (ModelTest-NG, IQ-TREE, DistanceCalculator). Assumes the normalised `alignment` and the imports from "Required Import and Normalisation" in `SKILL.md`.

## Distance Correction Models

For publication-grade pairwise distances, do NOT pick a model by rule of thumb. Run ModelTest-NG to select the best-fit substitution model by AIC/BIC, then apply that correction via IQ-TREE2 (`.mldist` output) or EMBOSS `distmat`. Hand-coded JC69 / K80 / blosum62 corrections via `Bio.Phylo.TreeConstruction.DistanceCalculator` are appropriate only for exploratory work.

```bash
modeltest-ng -i alignment.fasta -d nt -t ml
modeltest-ng -i alignment.fasta -d aa -t ml -p 4
```

Checked on ModelTest-NG 0.1.7: needs at least 3 real, non-ragged sequences (a 3-row toy example with `.`/lowercase characters, and A2M input with insert dots, both fail to parse — normalise and use a real alignment).

```python
from Bio.Phylo.TreeConstruction import DistanceCalculator

calculator = DistanceCalculator('blosum62')
distance_matrix = calculator.get_distance(alignment)
```

For substitution-model selection in the context of full phylogenetic inference, see `phylogenetics/modern-tree-inference`.
