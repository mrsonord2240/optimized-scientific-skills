## USRCAT (Ultra-Fast Shape Recognition + Atom Types)

USRCAT (Schreyer & Blundell 2012) extends Ultrafast Shape Recognition (USR) with atom-type information. Each molecule is represented as a 60-dimensional moment vector (12 moments × 5 atom types).

**Goal:** Encode a molecule into the 60-D USRCAT moment vector and score similarity against another molecule for alignment-free shape search.

**Approach:** Parse the SMILES, add hydrogens, generate one 3D conformer with ETKDGv3, compute RDKit USRCAT descriptors, and compare descriptor vectors with RDKit's USR score.

```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def usrcat(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    # 60 floats: 12 USR moments x 5 atom types
    # (all atoms, hydrophobic, aromatic, acceptor, donor)
    return rdMolDescriptors.GetUSRCAT(mol)

desc1, desc2 = usrcat('CCO'), usrcat('CCCO')
similarity = rdMolDescriptors.GetUSRScore(desc1, desc2)  # 1.0 for identical vectors
```

**Speed:** Descriptor calculation is linear in atoms and comparison is fixed-length, without pairwise alignment. Benchmark end-to-end throughput on the prepared conformer library before choosing a scale cutoff.

**Limit:** USRCAT is a coarse approximation. Predictive for analog identification; less precise for scaffold hopping.
