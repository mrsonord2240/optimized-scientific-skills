## Reactive Residue Taxonomy

| Residue | Nucleophile | Example compatible warheads | Design note |
|---------|-------------|-----------------------------|-------------|
| Cysteine | Thiol / thiolate | Acrylamide, haloacetamide, nitrile | Commonly targeted; local pKa and geometry control reactivity |
| Lysine | Amine | Sulfonyl fluoride, aldehyde | Aldehydes can form reversible imines with amines |
| Serine | Alcohol / alkoxide | β-lactam, boronate | Often requires catalytic activation |
| Threonine | Alcohol / alkoxide | Boronate | Context-dependent and less commonly targeted |
| Tyrosine | Phenol / phenolate | Sulfonyl fluoride, fluorosulfate | Local environment strongly affects reaction |
| Aspartate/Glutamate | Carboxylate | Residue-specific electrophiles require experimental validation | Do not infer aldehyde Schiff-base formation with carboxylates |

Cysteine is frequently targeted because its thiol/thiolate can be nucleophilic and its local environment can support selective proximity-driven reaction. GSH and off-target cysteines are competing thiols, not intrinsically distinguishable from the target by the warhead alone; the complete ligand's recognition, exposure, and intrinsic reactivity determine selectivity.

## Warhead Chemistry

| Warhead | SMARTS pattern | Reactivity | Reversibility | Cys-selective |
|---------|----------------|------------|---------------|----------------|
| Acrylamide | `[CX3](=[OX1])([NX3])[CX3]=[CX3]` | Moderate (Michael acceptor) | Usually irreversible | Often Cys-directed |
| Chloroacetamide | `[CX3](=[OX1])([NX3])[CH2]Cl` | High (SN2) | Irreversible | Often Cys-directed |
| Iodoacetamide (ABPP probes) | `[CX3](=[OX1])([NX3])[CH2][I]` | High (SN2) | Irreversible | Often Cys-directed |
| α-haloketone | `[CX3](=O)C[F,Cl,Br]` | Very high | Irreversible | Yes (but reactive) |
| Vinyl sulfone | `S(=O)(=O)C=C` | Moderate (Michael) | Irreversible | Yes |
| Sulfonyl fluoride | `S(=O)(=O)F` | Moderate | Irreversible | Lys/Tyr/Ser |
| Fluorosulfate (SuFEx) | `OS(=O)(=O)F` | Moderate | Irreversible | Tyr/Lys |
| Aldehyde | `[CX3H1](=O)` | Variable | Often reversible (covalent equilibrium) | Context-dependent Cys/Lys/Ser chemistry |
| Boronate (B-OH or B(OH)2) | `B(O)O` | Moderate | Reversible | Ser/Thr |
| Nitrile | `C#N` | Low | Reversible (Cys-S adduct) | Cys |
| Epoxide | `C1OC1` | High | Irreversible | Cys/Lys/Asp |
| α,β-unsaturated ketone | `[CX3](=O)C=C` | Moderate (Michael) | Irreversible | Cys |
| Isothiocyanate | `N=C=S` | High | Irreversible | Cys/Lys |
| Maleimide | `O=C1N(C(=O)C=C1)` | Often high | Commonly irreversible under assay conditions | Often Cys-directed |
| Cysteine-selective heterocycle | various | Moderate | Variable | Yes (designed) |

**Practical hierarchy:** Acrylamides are common attenuated electrophiles in cysteine-directed TCIs, including KRAS G12C, EGFR, and BTK programs. Haloacetamides are generally more intrinsically reactive, but actual selectivity must be measured for the complete molecule and target context.

**Overlapping SMARTS keys:** in `scripts/warhead_classifier.py`'s `WARHEAD_SMARTS`, `acrylamide`, `alpha_substituted_acrylamide`, and `methacrylamide` are not mutually exclusive by construction — an alpha-substituted or alpha-methyl acrylamide co-matches all three (verified on `C=C(C)C(=O)N1CCCCC1`). `classify_warheads()` correctly returns every matched key; a consumer that reduces this to one tier must pick the most specific matched key (`methacrylamide` > `alpha_substituted_acrylamide` > `acrylamide`), never a fixed lookup order, or it will silently report the least-specific tier.

## Decision Tree by Scenario

| Goal | Warhead choice | Reactivity tier |
|------|----------------|-----------------|
| Cysteine TCI program | Acrylamide is one common starting class | Measure complete-compound target and intrinsic reactivity |
| Cysteine probe program | Haloacetamides are one common class | Higher intrinsic reactivity can aid labeling but requires selectivity profiling |
| Lysine TCI (uncommon) | Sulfonyl fluoride | Moderate |
| Tyrosine TCI | Fluorosulfate (SuFEx) | Moderate |
| Reversible covalent | Warhead with demonstrated reversible adduct chemistry, such as selected cyanoacrylamides, aldehydes, nitriles, or boronates | Confirm reversibility experimentally |
| Activity-based protein profiling (ABPP) | Iodoacetamide / chloroacetamide | Very high |
| Boronic acid inhibitor (proteasome) | Boronate | Reversible |
| Aldehyde inhibitor (calpain) | Aldehyde | Reversible covalent |
