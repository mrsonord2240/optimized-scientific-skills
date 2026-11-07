## Ternary Complex Prediction Tools

| Tool | Approach | Strength | Fails when |
|------|----------|----------|------------|
| PRosettaC | Constrained PatchDock, RosettaDock refinement, PROTAC conformer generation, repacking, and clustering | PROTAC-specific published workflow | Performance varies by complex; Rosetta/service requirements |
| DeepTernary | Equivariant deep learning | Fast; SE(3) | OOD chemistry |
| AlphaFold3 | Unrestrained whole-complex prediction | Accepts proteins and ligands | No arbitrary user distance-restraint interface; benchmark PROTAC use |
| Boltz-1 / Boltz-2 | Unrestrained whole-complex prediction | Open local models | Limited PROTAC-specific validation |
| HADDOCK | Information-driven, restraint-guided docking | Mature integrative docking framework | Manual restraint specification |

**Decision:** Use a PROTAC-specific method such as **PRosettaC** for first-pass ternary modeling. AlphaFold3 or Boltz can provide unrestrained whole-complex predictions, but should be benchmarked on relevant ternary complexes. **DeepTernary** is released research code rather than a hosted API; validate it against relevant structures before prospective ranking.

**None of these tools run locally through this Skill.** All are external submissions or gated
local installs -- this Skill's only local, runnable content for ternary hypotheses is
`scripts/ternary_geometry_screen.py`, a linker-reach pre-filter (see "Ternary Complex Modeling
Workflow" below); it does not predict a structure or an interface score. To actually run one of
the methods above, read `references/ternary-prediction-tools.md` (per-tool URL, account, turnaround, output;
also a PRosettaC vs AlphaFold3 comparison).

Do not fabricate a numeric score from any of these tools when they have not actually been run --
state plainly that the step requires the external submission above.
