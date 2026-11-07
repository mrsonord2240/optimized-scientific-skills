# Illustrative BOIN and CRM designs
```bash
Rscript "$SKILL_DIR/scripts/dose.R" > dose-example.txt
```
The entry points are scripts/dose.R and scripts/map.R. Before running, set SKILL_DIR to the absolute skill directory and change into a newly created, empty working/output directory. Shell redirection truncates existing dose-example.txt or map-example.txt.

**Done when:** output includes the BOIN table and operating characteristics and CRM simulation under the script's stated dose curves, seeds and sample sizes. These are examples, not calibrated trial decisions; evaluate plausible scenarios and selection/overdose risk. Package APIs and versions must be checked.
For a single historical-data MAP example use `Rscript "$SKILL_DIR/scripts/map.R" > map-example.txt`. It is not an EXNEX or stratum-detachment model. See [dose and platform guidance](../references/platform-rar-dose-finding.md).
