## Docking Tool Taxonomy

| Tool | Scoring | Speed (sec/lig) | Best at | Fails when |
|------|---------|-----------------|---------|------------|
| AutoDock Vina 1.2 | Vina (empirical) | Hardware- and settings-dependent | Open, well-characterized baseline | Cross-dock; cryptic pockets; metal centers |
| SMINA | Vina + flexible + custom | Hardware- and settings-dependent | Custom scoring; flexible side chains | Same Vina-scoring caveats |
| Vinardo | Modified Vina scoring | Hardware- and settings-dependent | Alternative empirical score | Validate on target-relevant controls |
| GNINA 1.1 | CNN or Vina scoring | GPU- and settings-dependent | CNN-assisted pose ranking | Validate transfer to the target and chemotype |
| AutoDock 4 | AD4 + grid maps | Hardware- and settings-dependent | Legacy reference | More setup than Vina |
| DOCK 6/7 | DOCK + Amber | Hardware- and settings-dependent | UCSF DOCK ecosystem | Steep learning curve |
| Glide (Schrodinger) | GlideScore | License and hardware-dependent | Commercial docking workflow | License cost |
| GOLD (CCDC) | GOLDScore / ChemScore | License and hardware-dependent | Commercial workflow; metal options | License cost |
| FlexX (BioSolveIT) | FlexX | License and hardware-dependent | Fragment-based placement | License cost |
| rDock | rDock | Hardware- and settings-dependent | Open-source alternative | Validate maintenance and target fit |
| DiffDock-L | Diffusion-generative | GPU- and settings-dependent | Pose sampling for cross-docking | Validate geometry with PoseBusters; see ml-docking-rescoring |
| EquiBind | Equivariant NN | GPU- and settings-dependent | Single-shot pose generation | Requires independent geometry and ranking checks |
| Boltz-2 + GNINA rescore | Foundation model + CNN | GPU- and settings-dependent | Experimental multi-model workflow | Benchmark each evidence stream independently |

**Decision:** Use Vina as an open baseline and consider GNINA CNN rescoring when target-relevant redocking or cross-docking controls support it. For large libraries, calibrate a hierarchical Vina -> GNINA -> higher-cost follow-up workflow on measured enrichment, throughput, and retained chemotype diversity.

## Decision Tree by Scenario

| Scenario | Recommended workflow |
|----------|---------------------|
| Self-dock against known ligand pocket | GNINA `gnina --cnn_scoring rescore` (see `references/gnina.md`) |
| Cross-dock to apo or related-target structure | DiffDock-L pose + GNINA rescore + PoseBusters (failure modes: `references/failure-modes.md`) |
| Ultralarge library (10M+) | Calibrated hierarchical screen: property/alert triage -> Vina -> measured top fraction to GNINA -> higher-cost follow-up (see `references/ultralarge-screening.md`) |
| Cryptic pocket / induced fit | Receptor-ensemble docking and, where appropriate, a separately validated complex-prediction model |
| Allosteric / undefined site | P2Rank for pocket detection -> ensemble dock all pockets |
| Metal-coordinated ligand | GOLD (commercial) or manually parameterize Vina metal scoring |
| Covalent inhibitor | See `chemoinformatics/covalent-design`: DOCKovalent, HCovDock |
| Fragment screen (<300 Da) | rDock or constrained Vina with seed atoms |
| Hit-to-lead refinement | Use co-crystal structure if available; MD-relaxed receptor; FEP for affinity |
