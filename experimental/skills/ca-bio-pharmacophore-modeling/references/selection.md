## Pharmacophore Feature Types

| Feature | Common shorthand | Definition | Geometric tolerance |
|---------|------------|------------|----------------------|
| H-bond donor | D | -OH, -NH | 1.0-1.5 Å |
| H-bond acceptor | A | sp2 O / N (lone pair) | 1.0-1.5 Å |
| Hydrophobe | H | sp3 C / aromatic ring centroid | 1.5-2.0 Å |
| Aromatic ring | R | Aromatic ring centroid + normal | 1.0-1.5 Å |
| Positive ionizable | P | -NH3+, -NR3+ | 1.0-1.5 Å |
| Negative ionizable | N | -COO-, -SO3- | 1.0-1.5 Å |
| Halogen | X | Cl, Br, I (halogen bond donor) | 1.0-1.5 Å |
| Metal coordination | M | sp/sp2 N/O near metal | 0.5-1.0 Å |

Tolerances are pharmacophore-feature distance windows in the search. Tighter tolerances = fewer hits but more specific.

The ranges in this table are repository starting heuristics, not universal feature tolerances. Set final bounds from aligned-feature variability, coordinate uncertainty, and retrospective validation for the selected search engine.

The one-letter labels above are human-readable shorthand, not RDKit API codes. RDKit's shipped `BaseFeatures.fdef` uses family names such as `Donor`, `Acceptor`, `Hydrophobe`, `Aromatic`, `PosIonizable`, and `NegIonizable`. Its default feature definitions do not provide every halogen-bond or metal-coordination model; add and validate project-specific feature definitions when those interactions matter.

## Method Taxonomy

| Method | Origin | Use case | Fails when |
|--------|--------|----------|------------|
| Ligand-based (LBP) | Catalyst, MOE, RDKit Pharm3D | Multiple actives, no crystal | <3 actives; flexible actives |
| Receptor-based (RBP) | apo2ph4, LigandScout, PLIP | Co-crystal or a defined apo pocket | Uncertain pocket conformation |
| Common pharmacophore | Validated alignment/feature-consensus workflow; RDKit can represent and query the resulting model | Consensus from active set | Diverse actives or uncertain bioactive conformers confound alignment |
| Pocket-conditioned generation (PharmacoForge) | Flynn et al. 2025 | Generate candidate pharmacophores from a protein pocket | Does not directly generate molecules; pretrained model required |
| Active learning pharmacophore | Catalyst variant | Iterative refinement | Custom; not standard |

## Decision Tree by Scenario

| Scenario | Method | Tools |
|----------|--------|-------|
| Co-crystal structure available | Interaction-derived receptor model | PLIP or LigandScout + Pharmit |
| Apo structure with a defined pocket | Apo receptor model | apo2ph4; export LigandScout PML |
| Multiple active compounds, no crystal | Ligand-based common pharmacophore | Alignment plus consensus-feature derivation in validated custom or external tooling; RDKit Pharm3D can apply the resulting model |
| Single active compound | Single-conformer pharmacophore | RDKit Pharm3D from bioactive conformer |
| Scaffold hopping prospective | Receptor-based + shape filter | apo2ph4 or interaction-derived model + shape search |
| Cross-target SAR transfer | Common pharmacophore across targets | Manual + LigandScout |
| Generate pocket-conditioned pharmacophores | PharmacoForge | Diffusion model followed by library retrieval |
| Library pre-filtering | Pharmacophore screen | Pharmit search |
