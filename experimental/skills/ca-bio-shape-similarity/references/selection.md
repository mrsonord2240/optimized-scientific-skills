## Shape Method Taxonomy

| Tool | Speed | Approach | Open-source | Fails when |
|------|-------|----------|-------------|------------|
| ROCS / FastROCS (OpenEye) | Hardware/database/conformer-dependent; vendor reports millions of conformers/s for FastROCS | Gaussian shape + color | No | License and prepared database |
| ROCS X | Trillion-scale reaction/synthon space on Orion | FastROCS plus Bayesian-bandit sampling | No | Commercial cloud workflow |
| USRCAT | Very fast alignment-free descriptor comparison | Moment-based + atom types | Yes | Coarse approximation |
| Open3DAlign (RDKit) | Medium | MMFF atom-type/charge-weighted alignment | Yes | Requires compatible typed 3D structures |
| ShaEP | Benchmark on actual conformers/hardware | Field-based (shape + ESP) | Free binary; inspect license | Requires valid 3D structures and charges for ESP |
| ESPSim | Benchmark on actual workload | Electrostatic + shape | Yes | Limited public benchmarks |
| Phase-Shape (Schrödinger) | commercial | Shape + pharmacophore | No | Commercial |
| USR (original) | Very fast alignment-free comparison | Moment-based only | Yes | No atom-type information |

**Decision:** Select a shape method by matched retrieval/enrichment performance, conformer preparation, throughput, licensing, and score semantics. USRCAT is useful as a fast prefilter; Open3DAlign provides an open alignment method; ROCS/FastROCS provide commercial shape/color workflows.

## Decision Tree by Scenario

| Scenario | Method | Notes |
|----------|--------|-------|
| Large prepared library | USRCAT pre-filter + Open3DAlign rescore | Choose rescore budget from measured retrieval saturation; see `references/usrcat.md`, `references/open3dalign.md` |
| Production VS for scaffold hop | ROCS + color (commercial) | Industry standard |
| Scaffold hopping prospective | Open3DAlign with conformer ensemble | Shape + flexibility; see `references/open3dalign.md` |
| Bioisostere replacement | ROCS color with neutral scoring | Pharmacophore-equivalent matches |
| Patent space carve-out | Shape constraint + 2D dissimilarity | Combine shape + dissimilar scaffold |
| Library diversity assessment | USRCAT k-nearest neighbor | Fast; see `references/usrcat.md` |
| Crystal-bound conformer template | Open3DAlign starting from co-crystal pose | Bioactive shape; see `references/open3dalign.md` |
| Cross-target screening | Shape + pharmacophore feature | Combined screen |
| ESP-relevant pocket, electrostatic bioisosteres | ShaEP or ESPSim | See `references/esp-similarity.md` |

## Tanimoto-Combo Scoring (ROCS Standard)

TanimotoCombo = Tanimoto_shape + Tanimoto_color

- Tanimoto_shape: volume overlap normalized
- Tanimoto_color: pharmacophore feature overlap

Each component is normalized from 0 to 1, so TanimotoCombo ranges from 0 to 2. It is a sum, not an average. Select follow-up thresholds from a relevant benchmark or enrichment study; a single cutoff is not portable across query preparation, color-force-field settings, and library composition.
