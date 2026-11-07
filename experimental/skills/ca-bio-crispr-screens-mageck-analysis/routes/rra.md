# Compare two conditions
~~~bash
mageck test --count-table screen.count.txt --treatment-id Drug_r1,Drug_r2 --control-id Veh_r1,Veh_r2 --norm-method median --output-prefix drug_vs_veh --gene-lfc-method median --sort-criteria pos
~~~
Use matched vehicle for drug-specific effects. For paired samples, use --paired with treatment and control ordered by matching pair. For heavy selection use --norm-method control --control-sgrna ntcs.txt after reviewing controls. Inspect both neg and pos FDR and LFC, mapping and guide evidence. Source conventions include FDR <0.05 and absolute LFC >1; they are not proof of validation.
Done when contrast, normalization, guide-level results and both directions are reviewed.
