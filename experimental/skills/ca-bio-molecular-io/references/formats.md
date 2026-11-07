## Format Taxonomy

| Format | Dim | Stereo | Charges | Strength | Fails when |
|--------|-----|--------|---------|----------|------------|
| SMILES | 2D | Atom chirality `@/@@`; double-bond `/` and `\` | Atom-local formal charges | Compact, web-friendly, fast parse | Loses absolute coordinates; aromatic perception ambiguous across toolkits; tautomers not canonical |
| InChI | 2D | `/b`, `/t`, `/m`, `/s` stereo sublayers | `/q` charge and `/p` added/removed-proton sublayers; `/p` is not a pH model | Canonical by construction; cross-toolkit identity | Standard InChI normalizes mobile-H forms; limited organometallic stereo; large molecules may require special handling |
| SDF V2000 | 2D/3D | Wedge bonds | M CHG line | Industry default; metadata via tags | 999-atom count-field limit; cannot encode multi-component reactions; query atoms ambiguous |
| SDF V3000 | 2D/3D | Wedge + stereo flag | Inline charge | No atom limit; query support; rich properties | Some software (legacy) cannot read; verbose |
| MOL2 (Tripos) | 3D | Common records rely on 3D coordinates and toolkit perception; no portable explicit stereo field | Per-atom partial | SYBYL atom types preserved for docking | Atom-type dialects diverge (SYBYL vs Corina); RDKit MOL2 parser brittle |
| PDB | 3D | None | None standard | Universal protein format | No bond orders; aromatic perception lost; ligand names truncated to 3 chars |
| PDBQT | 3D | None | Gasteiger / AD4 | AutoDock-ready; torsion tree encoded | Specific to docking; no aromaticity layer |
| BinaryCIF (MMTF retired) | 3D | Encoded | Encoded | BinaryCIF (`.bcif`) is the current compact structural format; RCSB stopped serving MMTF files on July 2, 2024 and recommends BinaryCIF (RCSB PDB 2024) | Not all toolkits parse; binary format |
| CDX/CDXML | 2D | Drawing | Drawing | ChemDraw native | Not a structural format; converts unreliably |
| InChIKey | Hash | Stereo layer | n/a | Database key, fast lookup | Collision probability depends on key block and collection size; cannot recover structure |

## Aromaticity Perception (most common silent error)

Different toolkits perceive aromaticity differently. The same SMILES round-tripped between toolkits may produce different canonical strings and different fingerprints.

| Model | Toolkit | Rule | Symptom of mismatch |
|-------|---------|------|---------------------|
| Daylight | OpenEye, Daylight | 4n+2 π on planar ring | Furan, thiophene aromatic |
| RDKit default | RDKit | Daylight-like with extensions for fused / N-containing | Compatible with Daylight for drug-like molecules |
| MDL | Available in several toolkits, including RDKit as `AROMATICITY_MDL` | Five-membered rings are not aromatic unless part of a fused aromatic system; only C/N and one-electron donors qualify; exocyclic double bonds exclude an atom | Five-membered heteroaromatics and exocyclic-bond systems may differ from the default RDKit model |
| OpenEye | OEAroModel | Several modes | Charged thiophene non-aromatic in MDL but aromatic in OpenEye |

**Fix:** Always re-canonicalize via the toolkit doing analysis. If aromaticity must be reassigned explicitly in RDKit, use a concrete model, for example `Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_RDKIT)`, after the molecule is in an appropriate sanitized or kekulized state.

## Stereochemistry Layers

Stereo loss is the second most common silent error. Each format encodes stereo differently:

- SMILES: `@/@@` for tetrahedral, `/` and `\` for cis/trans double bonds
- InChI: `/b` for double-bond stereo, `/t` for tetrahedral stereo, and `/m` plus `/s` for inversion/overall stereo type
- SDF: wedge/hash bond + parity 0/1/2; cis/trans encoded via bond direction
- MOL2: common Tripos records provide 3D coordinates but no portable wedge or explicit stereo field; verify toolkit perception by round trip

**Round-trip tests:** If `Chem.MolToSmiles(Chem.MolFromSmiles(smi), isomericSmiles=True)` does not preserve the represented stereochemistry, inspect whether the source contained stereo markers and whether any step called `Chem.RemoveStereochemistry()` or discarded stereochemical coordinates/bond directions. Sanitization alone does not intentionally remove valid stereochemistry. If `MolFromMolFile` returns a molecule missing wedge bonds, inspect the source's coordinates, bond directions, and parity encoding.
