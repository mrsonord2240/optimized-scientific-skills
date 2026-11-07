## Pharmacophore Search (Pharmit / Pharmer)

For library screening, configure feature types, centers, radii, and optional shape constraints in Pharmit, or use a Pharmer database and query produced in the format required by the installed release. Do not pass LigandScout PML or a vendor `.ph4` file directly unless the selected interface documents that import path. Pharmit reported searching millions of conformers in seconds to minutes; actual runtime depends on query selectivity, database size, and deployment (Sunseri & Koes 2016).

## Pharmacophore Quality Validation

Evaluate a pharmacophore by:

1. **Retrospective enrichment**: a stated metric on target-relevant actives and inactives/decoys. DUD-E can provide a benchmark with known decoy-construction biases; COCONUT is a natural-products collection, not a target-specific active/decoy benchmark.
2. **Geometric tightness**: feature distance variance across actives
3. **Selectivity**: false positives in inactive set should be low
4. **Specific consistency**: pharmacophore matches each active's bioactive conformer

```python
def pharmacophore_enrichment(query_pharmacophore, actives, inactives,
                             matches_pharmacophore):
    """Return active/inactive match-rate enrichment for a supplied matcher."""
    if not actives or not inactives:
        raise ValueError('actives and inactives must both be non-empty')
    n_active_match = sum(
        bool(matches_pharmacophore(mol, query_pharmacophore))
        for mol in actives)
    n_inactive_match = sum(
        bool(matches_pharmacophore(mol, query_pharmacophore))
        for mol in inactives)
    active_rate = n_active_match / len(actives)
    inactive_rate = n_inactive_match / len(inactives)
    return float('inf') if inactive_rate == 0 else active_rate / inactive_rate
```

For this repository, enrichment >=5x may be used as a starting triage heuristic only after the active/decoy construction and matching policy are documented. Report the full metric and uncertainty, and calibrate the acceptance threshold on the project dataset.

## Pocket-Conditioned Pharmacophore Generation (PharmacoForge)

PharmacoForge (Flynn et al. 2025) applies a diffusion model to a protein pocket and generates candidate 3D pharmacophores. It does **not** directly generate molecular structures from an input pharmacophore. The validated workflow is:

1. Prepare the protein pocket in the representation required by the published PharmacoForge release.
2. Sample and rank pocket-conditioned pharmacophores.
3. Convert a selected pharmacophore into the query representation used by the search engine.
4. Retrieve matching, purchasable compounds and evaluate them with docking, strain, and physical-validity checks.

The paper compares pharmacophore and downstream retrieval performance with other pocket-based approaches; it does not support a drug-likeness or novelty comparison with REINVENT.
