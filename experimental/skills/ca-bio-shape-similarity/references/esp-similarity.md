## ESP Similarity (Electrostatic)

ShaEP and ESPSim extend shape with electrostatic surface potential overlap. For ESP-relevant pharmacophores (binding pockets with strong electrostatics):

### ShaEP (external binary, field-based)

RDKit has no Mol2 writer, so build the `.mol2` inputs ShaEP requires with Open Babel first (see `chemoinformatics/molecular-io` for other conversions):

```bash
obabel -:"CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1" -O query.mol2 --gen3D
obabel -:"CC(=O)Nc1ccc(C(=O)c2ccc(F)cc2)cc1" -O target.mol2 --gen3D
shaep -q query.mol2 target.mol2 -s aligned_hits.sdf similarity.txt
```

`similarity.txt` reports `shape_similarity`, `ESP_similarity`, and their average per target -- checked on ShaEP 1.4.2 (`-q` query, `-s` superimposed structures; with no `--output-file`, the last filename in the input list is the similarity output; ShaEP is a separate binary, not part of RDKit).

`obabel --gen3D` can print `NaN in calculated coordinates` on fused-ring scaffolds (seen on a naphthalene) and still write a usable mol2. Before trusting the ShaEP score, count zero and repeated coordinate rows in the mol2 `@<TRIPOS>ATOM` block; a good file gives `zero=0 dup=0`, a failed build gives `zero` equal to the atom count:

```bash
awk '/^@<TRIPOS>ATOM/{f=1;next} /^@<TRIPOS>/{f=0} f&&NF>=6{n++; k=$3" "$4" "$5; if(k in s)d++; s[k]=1; if($3+0==0&&$4+0==0&&$5+0==0)z++} END{print "atoms="n+0,"zero="z+0,"dup="d+0}' target.mol2
```

### ESPSim (RDKit-native, Python)

ESPSim (`pip install espsim`) embeds, aligns, and scores shape and ESP similarity without leaving RDKit:

```python
from rdkit import Chem
from espsim import EmbedAlignScore

query = Chem.AddHs(Chem.MolFromSmiles('CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1'))
target = Chem.AddHs(Chem.MolFromSmiles('CC(=O)Nc1ccc(C(=O)c2ccc(F)cc2)cc1'))

# prbNumConfs/refNumConfs: conformers generated internally for alignment search
shape_sim, esp_sim = EmbedAlignScore(target, [query], prbNumConfs=10, refNumConfs=10)
```

`shape_sim` is a shape Tanimoto in [0, 1]. `esp_sim` uses the default `metric='carbo'` (Carbo similarity index) and is **not** bounded to [0, 1] -- it can go negative for anti-correlated electrostatic potentials. Pass `renormalize=True` to rescale it to [0, 1] if a bounded score is needed. Default partial charges are RDKit Gasteiger; pass `prbCharge`/`refCharge` for higher-accuracy charges (e.g. from a QM calculation).

ESP scoring catches electrostatic-equivalent bioisosteres that pure shape misses (carboxylate vs tetrazole same charge).
