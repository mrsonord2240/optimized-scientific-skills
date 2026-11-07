## Binding Site Detection

When the binding pocket is not known (apo target, novel allosteric site):

| Tool | Approach | Output |
|------|----------|--------|
| P2Rank (Krivak 2018) | ML on protein surface descriptors | Ranked pocket list with center coords |
| fpocket (Le Guilloux 2009) | Voronoi tessellation | Pocket descriptor list |
| DoGSiteScorer | Geometric + drugability | Pocket list with score |
| AutoSite (Vina) | Affinity map clustering | Pocket centers |
| AlphaFill | Transplant ligands/cofactors from homologous experimental structures into AlphaFold models | Plausible binding-site components for review |

```bash
prank predict -f receptor.pdb -o pockets/
```

**On Windows, `prank.bat` returns immediately with no output when invoked through a shell
wrapper** (`cmd.exe /c`, or a bash script calling it) -- confirmed on this fixture: it printed
only the `cmd.exe` banner and produced no `pockets/` directory. Call the jar directly instead,
which does run and predict:

```bash
java -cp "<p2rank_dir>/bin/p2rank.jar;<p2rank_dir>/bin/lib/*" cz.siret.prank.program.Main \
     predict -f receptor.pdb -o pockets/
```

Verified on PDB 3PTB: 4 pockets in ~4 s, `pockets/rec.pdb_predictions.csv` with `pocket1` centered
at (-1.52, 14.47, 17.47) -- the same box center used throughout this Skill's Vina examples, which
sits on the real benzamidine site. P2Rank output `<receptor>_predictions.csv` lists pocket centers
with scores. The highest model score does not identify a pocket as orthosteric or biologically
relevant; verify ranked pockets against co-crystal, mutagenesis, SAR, or other structural evidence.
