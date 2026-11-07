## Receptor-Based Pharmacophore (apo2ph4 workflow)

**Goal:** Derive a pharmacophore from a protein binding-pocket structure without requiring a bound ligand.

**Approach:** Identify donor, acceptor, and hydrophobic hot spots from apo-pocket geometry, cluster them, and assemble candidate pharmacophores. Heider et al. describe apo2ph4 in *J. Chem. Inf. Model.* 63:101-110 (2023). Use the source release's documented scripts and environment rather than assuming a packaged `apo2ph4` command: the published workflow writes LigandScout PML output, not a generic `.ph4` file. Treat conversion to Pharmit, Pharmer, MOE, or Phase as a separate, explicitly validated step because pharmacophore formats are not interchangeable.

When a co-crystal ligand is available, **derive pharmacophore directly from the ligand binding pose**: each ligand feature in contact with a complementary protein residue is part of the pharmacophore.

```python
from plip.basic import config
from plip.structure.preparation import PDBComplex

mol_complex = PDBComplex()
mol_complex.load_pdb('complex.pdb')
mol_complex.analyze()

for site in mol_complex.interaction_sets.values():
    for interaction in site.all_itypes:
        # Objects are interaction-class-specific. Inspect the documented fields
        # for HydrophobicContact, HydrogenBond, PiStacking, SaltBridge, etc.;
        # there is no universal `.type` or `.ligatom.coords` interface.
        interaction_class = type(interaction).__name__
        print(interaction_class, interaction)
```

PLIP exposes typed interaction records with class-specific ligand/protein atoms and coordinates. Map those records to pharmacophore features explicitly and retain the interaction class and source atom identifiers.
