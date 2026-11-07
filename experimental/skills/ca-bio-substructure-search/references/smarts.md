## SMARTS Grammar Essentials

| Token | Meaning | Example |
|-------|---------|---------|
| `[#6]` | Atom by atomic number | `[#6]` carbon (any hybridization) |
| `c` | Lowercase = aromatic | `c1ccccc1` benzene aromatic |
| `C` | Uppercase = aliphatic only | `C(=O)O` carboxylic acid carbon |
| `[CX4]` | Atom + connection count X | `[CX4]` sp3 carbon (4 connections) |
| `[CX3]=O` | Carbonyl (CX3 = sp2 with 3 bonds) | matches ketone, aldehyde, ester C |
| `[#6;R]` | Atom in ring | `[#6;R]` ring carbon |
| `[#6;!R]` | Atom not in ring | `[#6;!R]` acyclic carbon |
| `[#6;r6]` | Atom in 6-membered ring | `[#6;r6]` six-ring carbon |
| `[a]` | Any aromatic atom | `[a]` |
| `[!#1]` | Anything except H | `[!#1]` heavy atom |
| `[N;H2]` | N with exactly 2 H; neighboring chemistry unconstrained | `[NH2]` also matches non-amine `NH2` environments unless context is added |
| `[N+]` | Positively charged N | `[N+](=O)[O-]` nitro |
| `[$(...)]` | Recursive SMARTS | `[$(c1ccccc1)]` aromatic 6-ring atom |
| `[c]([F,Cl,Br,I])` | OR within brackets | aryl halide |
| `~` | Any bond type | `c~c` any aromatic-aromatic bond |
| `@` | Ring bond | `c@c` requires the matched bond to be in a ring |
| `-` | Single bond explicit | `C-C` |
| `=` | Double bond | `C=O` |
| `:` | Aromatic bond explicit | |

## Common SMARTS Patterns

| Pattern | SMARTS | Notes |
|---------|--------|-------|
| Hydroxyl (alcohol + phenol) | `[OX2H]` | OX2H avoids matching O- in OH- |
| Phenol only | `[OX2H][c]` | OH attached to aromatic carbon |
| Aliphatic OH only | `[OX2H][CX4]` | OH attached to sp3 C |
| Carboxylic acid | `[CX3](=O)[OX2H1]` | C(=O)OH |
| Carboxylate | `[CX3](=O)[O-]` | C(=O)O- (deprotonated) |
| Ester | `[CX3](=[OX1])[OX2][#6]` | C(=O)O-R; supports alkyl and aryl esters |
| Amide | `[CX3](=[OX1])[NX3]` | C(=O)N-R |
| Primary amine attached to carbon (excluding common amide-like N) | `[NX3;H2;$(N-[#6]);!$(N-[C,S,P]=[O,S,N])]` | Carbon-substituted -NH2; extend the exclusions for a project-specific amine definition |
| Secondary amine attached to two carbons | `[NX3;H1;$(N(-[#6])-[#6]);!$(N-[C,S,P]=[O,S,N])]` | Carbon-substituted -NH- excluding common amide-like N |
| Neutral tertiary amine attached to three carbons | `[NX3;H0;+0;$(N(-[#6])(-[#6])-[#6]);!$(N-[C,S,P]=[O,S,N])]` | Carbon-substituted -NR2 excluding common amide-like N |
| Quaternary amine | `[NX4+]` | -NR4+ |
| Nitro | `[N+](=O)[O-]` | -NO2 |
| Nitrile | `[CX2]#[NX1]` | -C#N |
| Sulfonamide | `[SX4](=[OX1])(=[OX1])[NX3]` | -S(=O)(=O)N |
| Aryl halide | `[c][F,Cl,Br,I]` | halogen on aromatic |
| Aliphatic halide | `[CX4][F,Cl,Br,I]` | halogen on sp3 C |
| Hydrogen bond donor | Use a named feature definition such as RDKit `BaseFeatures.fdef` or `Lipinski.NumHDonors` | `[#7,#8;!H0]` is only a simplified N/O-H query and is not a universal HBD model |
| Hydrogen bond acceptor | Use a named feature definition such as RDKit `BaseFeatures.fdef` or `Lipinski.NumHAcceptors` | No short universal SMARTS correctly captures every accepted HBA chemistry model |
| Michael acceptor | `[CX3]=[CX3][CX3]=O` | enone, acrylamide warhead |
| Aldehyde | `[CX3H1](=O)` | -CHO |
| Ketone | `[CX3;H0](=[OX1])([#6])[#6]` | Carbonyl carbon has two carbon substituents and no hydrogen |

## Basic Substructure Match

**Goal:** Test whether a molecule contains a SMARTS pattern and enumerate the matching atom indices.

**Approach:** Parse the molecule with `MolFromSmiles` and the pattern with `MolFromSmarts`, gate with `HasSubstructMatch`, then call `GetSubstructMatches` and map each atom index back to the molecule for inspection.

```python
from rdkit import Chem

def compile_smarts(smarts):
    if not isinstance(smarts, str) or not smarts.strip():
        raise ValueError('SMARTS must be a non-empty string')
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None or pattern.GetNumAtoms() == 0:
        raise ValueError(f'Invalid SMARTS: {smarts}')
    return pattern

mol = Chem.MolFromSmiles('c1ccc(O)cc1CCO')
pattern = compile_smarts('[OX2H]')

if mol.HasSubstructMatch(pattern):
    matches = mol.GetSubstructMatches(pattern)
    for match in matches:
        atoms = [mol.GetAtomWithIdx(i).GetSymbol() for i in match]
```

`HasSubstructMatch` returns bool, `GetSubstructMatches` returns tuple of tuples of atom indices.

## Recursive SMARTS for Context-Aware Patterns

`[$(pattern)]` matches an atom that *also* matches the entire pattern starting from itself. Critical for context-aware matching.

```python
# Aromatic carbon attached to a carbonyl
pat = Chem.MolFromSmarts('[$(c[C](=O))]')

# Aniline-type N (aromatic carbon-N-H)
pat = Chem.MolFromSmarts('[$([NX3;H2][c])]')

# Neutral tertiary amine with three sp3-carbon neighbors
pat = Chem.MolFromSmarts('[$([NX3]([CX4])([CX4])[CX4])]')

# H-bond donor (per Lipinski, exclude quaternary)
hbd = Chem.MolFromSmarts('[#7,#8;!H0;!$([NX3+])]')

# For H-bond acceptors, use RDKit's maintained Lipinski/feature definitions
# instead of an ad hoc universal SMARTS.
from rdkit.Chem import Lipinski
n_acceptors = Lipinski.NumHAcceptors(mol)
```
