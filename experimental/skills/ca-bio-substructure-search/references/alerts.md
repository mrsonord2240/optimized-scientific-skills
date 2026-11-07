## Structural-Alert Filter Catalogs

| Filter | Origin | Patterns | Use case | Failure mode |
|--------|--------|----------|----------|--------------|
| PAINS_A | Baell & Holloway 2010 | 16 | Most populated source-data patterns (>=150 analogues per pattern) | Many false positives in primary screens; legitimate medicines flagged |
| PAINS_B | Baell & Holloway 2010 | 55 | Intermediate source-data population (15-149 analogues per pattern) | Similar |
| PAINS_C | Baell & Holloway 2010 | 409 | Least populated source-data patterns (1-14 analogues per pattern) | Most permissive |
| BRENK | Brenk 2008 (DDS unsuitable) | 105 | Reactive / toxicity / undesirable | Useful for fragment / virtual library |
| NIH | NIH MLSMR | 180 in RDKit 2024.09 | Reactive groups, unstable | Legacy filter; verify count after toolkit upgrades |
| ZINC | ZINC clean-leads | 50 in RDKit 2024.09 | Drug-like cleanup | Verify definitions after toolkit upgrades |
| Glaxo / Eli Lilly | Vendor lists | varies | Internal "ugly" filters | Often unpublished |
| REOS | Walters & Murcko 2002 | property + structural | Requires an explicit project implementation; not bundled here | Hand-curated thresholds |

The PAINS A/B/C families encode pattern population in the original screening dataset, not increasing or decreasing external evidence strength.

## When to Apply Each Filter

| Scenario | Catalog | Reason |
|----------|---------|--------|
| Hit validation from biochemical screen | PAINS_A | Identify assay-interference candidates |
| Library prep for HTS | PAINS_A + Brenk + ZINC | Remove clearly bad |
| Fragment library design | Brenk + ZINC | Remove reactive; PAINS less critical at fragments |
| Lead optimization | None mandatory | Filters can exclude valid leads |
| Natural product analog | None | Filters trained on synthetic chemistry |
| Covalent inhibitor design | Skip warhead filter | Warheads ARE the design |

**Critical:** Capuzzi et al. (2017) found PAINS alerts in 87 FDA-approved small-molecule drugs. PAINS is a *flag for assay validation*, not a *killing filter*. Catalog rates are collection-dependent: in the reference 800-compound lead-like collection, BRENK flagged 36% and PAINS_A+BRENK+ZINC flagged 39.4% in union. Review BRENK hits by alert category rather than treating a high flag rate as a deletion quota.

## PAINS Filter

**Goal:** Split a molecule list into PAINS-flagged and PAINS-clean sets using one or more PAINS catalog tiers.

**Approach:** Configure `FilterCatalogParams` with the requested catalog enums, build a `FilterCatalog` once, and for each molecule use `GetFirstMatch` to either bucket it as clean or record the matching pattern description.

```python
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

def pains_filter(mols, catalogs=('PAINS_A',)):
    params = FilterCatalogParams()
    for cat in catalogs:
        params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, cat))
    catalog = FilterCatalog(params)

    flagged = []
    clean = []
    invalid = []
    for index, mol in enumerate(mols):
        if mol is None:
            invalid.append(index)
            continue
        entry = catalog.GetFirstMatch(mol)
        if entry is None:
            clean.append(mol)
        else:
            flagged.append((mol, entry.GetDescription()))
    return clean, flagged, invalid
```

Available catalog names: `PAINS_A`, `PAINS_B`, `PAINS_C`, `PAINS` (all), `BRENK`, `NIH`, `ZINC`, `ALL`. `invalid` contains zero-based input positions whose molecules did not parse; report them rather than silently dropping records.

## Reaction-Reactive Group Filter (custom)

For HTS triage, filter electrophilic warheads (acrylamide, chloroacetamide, etc.) unless designing covalent inhibitors.

**Goal:** Flag molecules containing electrophilic warheads or other reactive functional groups that would interfere with biochemical HTS.

**Approach:** Maintain a named SMARTS dictionary of reactive groups (acid halides, epoxides, Michael acceptors, etc.), then per molecule scan each pattern with `HasSubstructMatch` and return the first matching warhead name.

```python
REACTIVE_SMARTS = {
    'acid_anhydride': '[CX3](=O)O[CX3](=O)',
    'acid_halide': '[CX3](=O)[F,Cl,Br,I]',
    'alpha_halo_carbonyl': '[CX3](=O)C([F,Cl,Br,I])',
    'aldehyde_reactive': '[CX3H1](=O)[#6;X4]',  # aliphatic aldehydes
    'epoxide': 'C1OC1',
    'aziridine': 'C1NC1',
    'isocyanate': '[NX2]=C=[OX1]',
    'isothiocyanate': '[NX2]=C=[SX1]',
    'beta_lactam': 'C1(=O)NCC1',
    'sulfonyl_halide': '[SX4](=O)(=O)[F,Cl,Br,I]',
    'Michael_acceptor': '[CX3]=[CX3][CX3]=O',
    'vinyl_sulfone': '[SX4](=O)(=O)C=C',
}

# Compile once per filter configuration, not once per molecule.
COMPILED_REACTIVE_SMARTS = {
    name: compile_smarts(smarts) for name, smarts in REACTIVE_SMARTS.items()
}

def reactive_filter(mol, exclude_warheads=True):
    if not exclude_warheads:
        return False, None
    if mol is None:
        return False, None
    for name, pattern in COMPILED_REACTIVE_SMARTS.items():
        if mol.HasSubstructMatch(pattern):
            return True, name
    return False, None
```

For covalent-inhibitor design, see `chemoinformatics/covalent-design`; these warheads are the desired chemistry, not noise to filter. Treat custom reactive alerts as flags, not deletions: broad Michael-acceptor patterns also flag designed enones (including curcumin), and beta-lactam patterns flag antibiotics. Review assay format, intended mechanism, and orthogonal evidence before excluding a series.
