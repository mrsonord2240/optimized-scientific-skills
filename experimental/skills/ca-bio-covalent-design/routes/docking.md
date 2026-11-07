# Covalent docking

```bash
autogrid4 -p target.gpf -l target.glg && autodock4 -p ligand_target.dpf -l ligand_target.dlg
```

This command assumes prepared AutoDock4 covalent flexible-receptor inputs; follow `references/docking.md` for the source workflow and tool constraints. Inspect the reaction-specific electrophile and nucleophile geometry, and do not infer labeling from docking alone. Failure modes are in `references/troubleshooting.md`.

Done when the requested prepared-input docking result and evidence limits are recorded.
