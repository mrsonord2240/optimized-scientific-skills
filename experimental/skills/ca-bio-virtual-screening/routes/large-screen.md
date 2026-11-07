# GNINA and large screens

```bash
gnina -r receptor.pdb -l ligand.sdf --autobox_ligand reference_ligand.sdf --cnn_scoring rescore -o poses.sdf.gz --num_modes 9 --exhaustiveness 8
```

Confirm installed GNINA modes and validate scoring on target-relevant controls; see `references/gnina.md`. The source ultralarge orchestrator is a skeleton with unimplemented project helpers, not a runnable pipeline; see `references/ultralarge-screening.md` for design and throughput caveats.

Done when actual scoring outputs and validation are recorded.
